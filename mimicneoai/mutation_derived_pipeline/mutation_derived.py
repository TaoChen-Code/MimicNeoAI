# coding=utf-8
"""
Launcher for the mutation-derived (neoantigen) pipeline:
- QC (fastp)
- Variant calling (WES/WGS/RNA where supported)
- Annotation (VEP)
- HLA typing
- Neoantigen identification via pVACseq
"""

from __future__ import annotations

import argparse
import multiprocessing.pool
import os
import traceback
from multiprocessing import Manager
from pathlib import Path
from typing import Any, Dict, List, Optional

from mimicneoai.functions.fastp import fastp
from mimicneoai.functions.hlatyping import hlahd
from mimicneoai.functions.pipline_tools import tools
from mimicneoai.mutation_derived_pipeline.scripts.annotation import annotation_vcf
from mimicneoai.mutation_derived_pipeline.scripts.hla_binding_pred import Pvacseq
from mimicneoai.mutation_derived_pipeline.scripts.variants_calling import (
    variants_calling_start,
)

# -------- Constants --------
FLAG = "Neoantigen"
STEP_NAME = {
    "QC": "00.QC",
    "alignment": "01.alignment",
    "sort": "02.sort",
    "rmdup": "03.rmdup",
    "bqsr": "04.BQSR",
    "variants_calling": "05.variants_calling",
    "vqsr": "06.vqsr",
    "annotation": "07.vep",
    "hla": "08.hlatyping",
    "pvacseq": "09.binding_prediction",
}


# -------- Non-daemon Pool (only needed if workers spawn child processes) --------
class _NoDaemonProcess(multiprocessing.Process):
    """Process class that always behaves as non-daemon (allows child processes)."""

    # For Python < 3.8
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


if hasattr(multiprocessing, "get_start_method"):
    # Python ≥ 3.8: override Pool's Process to use the non-daemon variant
    class NoDaemonPool(multiprocessing.pool.Pool):
        @staticmethod
        def Process(_, *args, **kwargs):
            return _NoDaemonProcess(*args, **kwargs)
else:
    # Python < 3.8 fallback
    class NoDaemonPool(multiprocessing.pool.Pool):
        Process = _NoDaemonProcess


# -------- Worker helpers --------
def _variants_calling_and_annotation(
    sample: str,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """Run variant calling and optional annotation for a single sample."""
    host_variants_calling = configure["others"]["host_variants_calling"]
    do_annotation = configure["others"]["annotation"]
    species = configure["others"]["species"]  # e.g., "human"
    seq_type = configure["others"]["seq_type"]  # e.g., "wes"
    tumor_with_matched_normal = configure["others"]["tumor_with_matched_normal"]

    # Dispatch to the appropriate variant-calling entry point
    if host_variants_calling:
        funcs = {
            "wes_human": variants_calling_start,
            # Extend here as more modalities/species are supported:
            # "wgs_human": wgs_human_mutation_calling_start,
            # "rna_human": rna_mutation_calling_start,
        }
        key = f"{seq_type}_{species}"
        func = funcs.get(key)
        if func is None:
            tool.write_log(f"No variant-calling function for key='{key}'", "error")
        else:
            func(sample, tool, configure, paths)

    # VEP-based annotation
    if do_annotation:
        if tumor_with_matched_normal:
            tumor_sample = sample.split(",")[0]
            _ = annotation_vcf(sample, tumor_sample, tool, paths, configure)
        else:
            _ = annotation_vcf(sample, sample, tool, paths, configure)


def _start_one_sample(
    sample: str,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """End-to-end per-sample workflow: QC → variants/annotation → HLA → pVACseq."""
    try:
        sample = str(sample)

        # Runtime toggles
        do_qc = configure["others"]["QC"]
        do_hlatyping = configure["others"]["hlatyping"]
        do_pvacseq = configure["others"]["peptides_identification_and_binding_prediction"]
        tumor_with_matched_normal = configure["others"]["tumor_with_matched_normal"]

        # 1) Read QC (fastp)
        if do_qc:
            if tumor_with_matched_normal:
                tumor_sample, normal_sample = sample.split(",")[0], sample.split(",")[1]
                fastp(sample, tumor_sample, configure, paths, tool)
                fastp(sample, normal_sample, configure, paths, tool)
            else:
                fastp(sample, sample, configure, paths, tool)

        # 2) Variant calling + annotation
        _variants_calling_and_annotation(sample, configure, paths, tool)

        # 3) HLA typing
        if do_hlatyping:
            if tumor_with_matched_normal:
                tumor_sample = sample.split(",")[0]
                hlahd(sample, tumor_sample, configure, paths, tool)
            else:
                hlahd(sample, sample, configure, paths, tool)

        # 4) Neoantigen prediction (pVACseq)
        if do_pvacseq:
            output_dir = configure["path"]["output_dir"]
            step_name_vep = configure["step_name"]["annotation"]
            step_name_hla = configure["step_name"]["hla"]
            pvacseq_runner = Pvacseq(tool)

            if tumor_with_matched_normal:
                tumor_sample = sample.split(",")[0]
                output_vep = f"{output_dir}/{tumor_sample}/{step_name_vep}/"
                output_hla = f"{output_dir}/{tumor_sample}/{step_name_hla}/"
                pvacseq_runner.run_pvacseq_parallel(
                    sample, tumor_sample, output_vep, output_hla, configure, paths
                )
            else:
                output_vep = f"{output_dir}/{sample}/{step_name_vep}/"
                output_hla = f"{output_dir}/{sample}/{step_name_hla}/"
                pvacseq_runner.run_pvacseq_parallel(
                    sample, sample, output_vep, output_hla, configure, paths
                )

    except Exception:
        # Keep error messages generic to avoid leaking sensitive environment details
        tool.write_log(f"Worker crashed:\n{traceback.format_exc()}", "error")


def _run_pipeline(
    samples: List[str],
    pool_size: int,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """Execute the pipeline across samples using a non-daemon process pool."""
    with NoDaemonPool(processes=pool_size) as pool:
        for sample in samples:
            pool.apply_async(
                _start_one_sample,
                (sample, configure, paths, tool),
                error_callback=tool.print_pool_error,
            )
        pool.close()
        pool.join()


def _peek_output_dir(cfg_path: str) -> Optional[str]:
    """Return configure.path.output_dir from YAML without fully parsing user config."""
    try:
        import yaml
    except Exception:
        return None
    try:
        with open(cfg_path, "r") as f:
            data = yaml.safe_load(f) or {}
        return (data.get("path", {}) or {}).get("output_dir", None)
    except Exception:
        return None


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point for the mutation-derived (neoantigen) pipeline."""
    parser = argparse.ArgumentParser(description="MimicNeoAI (Mutation-derived pipeline)")
    parser.add_argument(
        "-c", "--configure", required=True, help="Path to configuration YAML"
    )
    parser.add_argument(
        "-p", "--paths", dest="paths", required=False, help="Path to paths YAML"
    )
    parser.add_argument(
        "--workdir",
        default=None,
        help="Working directory root (default: config.path.output_dir or current directory)",
    )
    args = parser.parse_args(argv)

    # 1) Determine working directory precedence: --workdir > configure.path.output_dir > CWD
    cfg_output_dir = _peek_output_dir(args.configure)
    base_out = Path(cfg_output_dir)
    cfg_output_dir = str(base_out / "Mutation-derived")

    workdir = args.workdir or cfg_output_dir or os.getcwd()
    workdir = str(Path(workdir).resolve())

    # Initialize shared manager and logger/tool (log root = workdir)
    mgr = Manager()
    log_lock = mgr.Lock()
    tool_obj = tools(workdir, FLAG, log_lock)

    tool_obj.write_log(f"work_path: {tool_obj.sys_path}", "info")
    tool_obj.write_log(f"start_log: {tool_obj.start_log}", "info")
    tool_obj.write_log(f"cmd_log: {tool_obj.cmd_log_dir}", "info")
    if args.workdir is None and cfg_output_dir:
        tool_obj.write_log(
            f"workdir inferred from configure.path.output_dir: {cfg_output_dir}", "info"
        )

    # Load runtime configurations
    configure = tool_obj.get_configure(args.configure)
    paths = tool_obj.get_paths(args.paths)
    tool_obj.write_log(f"configures: {configure}", "info")
    tool_obj.write_log(f"paths: {paths}", "info")

    # Sample list and step-name mapping
    samples = [str(s) for s in configure["samples"]]
    configure = dict(configure)  # shallow copy to avoid mutating the original
    configure.setdefault("step_name", {}).update(STEP_NAME)
    base_out = Path(configure["path"]["output_dir"])
    configure["path"]["output_dir"] = str(base_out / "Mutation-derived")

    # Share variables across processes
    tool_obj.sharing_variable(mgr, samples)

    pool_size = int(configure["args"]["pool_size"])
    _run_pipeline(samples, pool_size, configure, paths, tool_obj)

    tool_obj.summary()
    return 0


if __name__ == "__main__":
    # Optional Windows support:
    # from multiprocessing import freeze_support; freeze_support()
    raise SystemExit(main())
