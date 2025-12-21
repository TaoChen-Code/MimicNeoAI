# coding=utf-8
from __future__ import annotations
import sys
import argparse
import os
import traceback
from multiprocessing import Manager
import multiprocessing.pool
from pathlib import Path
from typing import Any, Dict, List, Optional
import yaml

from mimicneoai.functions.fastp import fastp
from mimicneoai.functions.pipline_tools import tools
from mimicneoai.functions.hlatyping import hlahd
from mimicneoai.microbial_pipeline.scripts.microbial_peptides import (
    HostSequencesRemoving,
    VectorContaminationRemoving,
    MicrobialTaxasQuantification,
    MicrobialPeptidesIdentification,
    MicrobialPeptidesBindingPrediction,
)

FLAG = "MicrobialAntigen"
STEP_NAME = {
    "QC": "00.QC",
    "hg38": "01.HostSequencesRemovingStep1",
    "t2t": "02.HostSequencesRemovingStep2",
    "vector": "03.VectorContaminationRemoving",
    "pathseq": "04.MicrobialTaxaQuantificationStep1",
    "nucleic": "05.MicrobialTaxaQuantificationStep2",
    "blastx": "06.MicrobialPeptidesIdentification",
    "hla": "07.HlaTyping",
    "pvacbind": "08.MicrobialPeptidesBindingPrediction",
}

# -------- Non-daemon Pool (only needed if workers spawn child processes) --------
class _NoDaemonProcess(multiprocessing.Process):
    """Process whose daemon attribute is always False (allows nested pools)."""
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

if sys.version_info >= (3, 8):
    # Python 3.8+
    class NoDaemonPool(multiprocessing.pool.Pool):
        @staticmethod
        def Process(_, *args, **kwargs):
            # 丢掉第一个 ctx 参数，用自定义 _NoDaemonProcess
            return _NoDaemonProcess(*args, **kwargs)
else:
    # Python 3.7 early
    class NoDaemonPool(multiprocessing.pool.Pool):
        Process = _NoDaemonProcess


def _start_one_sample(
    sample: str,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """Run the microbial pipeline for a single sample with per-step guards."""
    o = configure["others"]

    run_qc                   = o["QC"]
    run_host_depletion       = o["run_host_depletion"]
    run_vector_decontam      = o["run_vector_decontamination"]
    run_pathseq              = o["run_pathseq"]
    run_microbial_peptide_identification    = o["run_microbial_peptide_identification"]
    run_hla_typing           = o["run_hla_typing"]
    run_binding_prediction   = o["run_binding_prediction"]

    if run_qc:
        try:
            fastp(sample, sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"QC error:\n{traceback.format_exc()}", "error")

    if run_host_depletion:
        try:
            HostSequencesRemoving(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"HostSequencesRemoving error:\n{traceback.format_exc()}", "error")

    if run_vector_decontam:
        try:
            VectorContaminationRemoving(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"VectorContaminationRemoving error:\n{traceback.format_exc()}", "error")

    if run_pathseq:
        try:
            MicrobialTaxasQuantification(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"MicrobialTaxasQuantification error:\n{traceback.format_exc()}", "error")

    if run_microbial_peptide_identification:
        try:
            MicrobialPeptidesIdentification(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"MicrobialPeptidesIdentification error:\n{traceback.format_exc()}", "error")

    if run_hla_typing:
        try:
            hlahd(sample, sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"HLA typing error:\n{traceback.format_exc()}", "error")

    if run_binding_prediction:
        try:
            MicrobialPeptidesBindingPrediction(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"Binding prediction error:\n{traceback.format_exc()}", "error")



def _run_pipeline(samples: List[str], pool_size: int, configure, paths, tool: tools) -> None:
    """Run the pipeline across samples using a (non-daemon) process pool."""
    with NoDaemonPool(pool_size) as pool:
        for sample in samples:
            pool.apply_async(
                _start_one_sample,
                (str(sample), configure, paths, tool),
                error_callback=tool.print_pool_error,
            )
        pool.close()
        pool.join()


def _peek_output_dir(cfg_path: str) -> Optional[str]:
    """Quickly read YAML to fetch path.output_dir; return None if missing or on error."""
    try:
        with open(cfg_path, "r") as f:
            cfg = yaml.safe_load(f) or {}
        return cfg.get("path", {}).get("output_dir")
    except Exception:
        return None


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point for the microbial antigen pipeline."""
    parser = argparse.ArgumentParser(description="Microbial Antigen pipeline")
    parser.add_argument("-c", "--configure", required=True, help="Path to configuration YAML")
    parser.add_argument("-p", "--paths", dest="paths", required=False, help="Path to paths YAML")
    parser.add_argument(
        "--workdir", default=None,
        help="Working directory (default: config.path.output_dir if set, else CWD)"
    )
    args = parser.parse_args(argv)

    # 1) Decide working/log root directory with precedence:
    #    --workdir (explicit) > configure.path.output_dir > current working directory
    cfg_output_dir = _peek_output_dir(args.configure)
    base_out = Path(cfg_output_dir)
    cfg_output_dir = str(base_out / "Microbial")
    workdir = args.workdir or cfg_output_dir or os.getcwd()
    workdir = str(Path(workdir).resolve())

    # 2) Initialize logger/tools (log root = workdir)
    mgr = Manager()
    log_lock = mgr.Lock()
    tool = tools(workdir, FLAG, log_lock)

    tool.write_log(f"work_path: {tool.sys_path}", "info")
    tool.write_log(f"start_log: {tool.start_log}", "info")
    tool.write_log(f"cmd_log: {tool.cmd_log_dir}", "info")
    if args.workdir is None and cfg_output_dir:
        tool.write_log(f"workdir inferred from configure.path.output_dir: {cfg_output_dir}", "info")

    # 3) Load configuration and paths (tools.get_configure also persists a copy for reproducibility)
    configure = tool.get_configure(args.configure)
    paths = tool.get_paths(args.paths)
    tool.write_log(f"configures: {configure}", "info")
    tool.write_log(f"paths: {paths}", "info")

    # 4) Fill in step names and collect sample list
    samples = list(configure.get("samples", []))
    configure = dict(configure)  # avoid mutating the original object
    configure.setdefault("step_name", {}).update(STEP_NAME)
    base_out = Path(configure["path"]["output_dir"])
    configure["path"]["output_dir"] = str(base_out / "Microbial")

    # 5) Initialize shared state & execute in parallel
    tool.sharing_variable(mgr, samples)
    pool_size = int(configure["args"]["pool_size"])
    _run_pipeline(samples, pool_size, configure, paths, tool)

    # 6) Final summary
    tool.summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
