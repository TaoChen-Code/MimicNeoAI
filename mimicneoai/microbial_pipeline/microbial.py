# coding=utf-8
from __future__ import annotations

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
    MicrobialTaxasQuantification,
    MicrobialPeptidesIdentification,
    MicrobialPeptidesBindingPrediction,
)

FLAG = "MicrobialAntigen"
STEP_NAME = {
    "QC": "00.QC",
    "hg38": "01.HostSequencesRemovingStep1",
    "t2t": "02.HostSequencesRemovingStep2",
    "pathseq": "03.MicrobialTaxasQuantificationStep1",
    "nucleic": "04.MicrobialTaxasQuantificationStep2",
    "blastx": "05.MicrobialPeptidesIdentification",
    "hla": "06.HlaTyping",
    "pvacbind": "07.MicrobialPeptidesBindingPrediction",
}

# -------- Non-daemon Pool (only needed if workers spawn child processes) --------
class _NoDaemonProcess(multiprocessing.Process):
    """Process whose daemon attribute is always False (allows nested pools)."""
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

if hasattr(multiprocessing, "get_start_method"):
    # Python ≥ 3.8: override Pool.Process to use _NoDaemonProcess
    class NoDaemonPool(multiprocessing.pool.Pool):
        @staticmethod
        def Process(_, *args, **kwargs):
            return _NoDaemonProcess(*args, **kwargs)
else:
    # Python < 3.8 fallback
    class NoDaemonPool(multiprocessing.pool.Pool):
        Process = _NoDaemonProcess


def _start_one_sample(
    sample: str,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """Run the microbial pipeline for a single sample with per-step guards."""
    QC = configure["others"]["QC"]
    host_sequences_removing = configure["others"]["host_sequences_removing"]
    microbial_taxas_quantification = configure["others"]["microbial_taxas_quantification"]
    microbial_peptides_identification = configure["others"]["microbial_peptides_identification"]
    hlatyping = configure["others"]["hlatyping"]
    microbial_peptides_bindingPrediction = configure["others"]["microbial_peptides_bindingPrediction"]

    if QC:
        try:
            fastp(sample, sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"QC error:\n{traceback.format_exc()}", "error")

    if host_sequences_removing:
        try:
            HostSequencesRemoving(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"HostSequencesRemoving error:\n{traceback.format_exc()}", "error")

    if microbial_taxas_quantification:
        try:
            MicrobialTaxasQuantification(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"MicrobialTaxasQuantification error:\n{traceback.format_exc()}", "error")

    if microbial_peptides_identification:
        try:
            MicrobialPeptidesIdentification(sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"MicrobialPeptidesIdentification error:\n{traceback.format_exc()}", "error")

    if hlatyping:
        try:
            hlahd(sample, sample, configure, paths, tool)
        except Exception:
            tool.write_log(f"HLA typing error:\n{traceback.format_exc()}", "error")

    if microbial_peptides_bindingPrediction:
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

    # 5) Initialize shared state & execute in parallel
    tool.sharing_variable(mgr, samples)
    pool_size = int(configure["args"]["pool_size"])
    _run_pipeline(samples, pool_size, configure, paths, tool)

    # 6) Final summary
    tool.summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
