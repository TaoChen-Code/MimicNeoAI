# coding=utf-8
"""
Launcher for the cryptic (lncRNA/sORF-derived) antigen pipeline.
It mirrors the original run.sh flow:

00  QC (tumor + optional control)
01  STAR alignment
02  lncRNA/sORF discovery (known + novel/trinity)
03  Salmon quant (tumor)
04  Salmon quant (control)
05  HLA typing (hlahd)
06  Extract aberrantly expressed sORF peptides (aeSEPs)
07  HLA binding prediction (pvacbind/IEDB)

This module:
- Uses absolute imports within the 'mimicneoai' package
- Avoids side effects at import time
- Reads YAML 'configure' and 'paths' like other pipelines
- Spawns a process pool across samples
"""

from __future__ import annotations
import os
from pathlib import Path
import yaml
import argparse
import traceback
import shlex
import sys
from typing import Dict, Any, List, Tuple
from multiprocessing import Manager
import multiprocessing.pool

from importlib.resources import files

from mimicneoai.functions.pipline_tools import tools


# ---------------------- Constants ----------------------
FLAG = "CrypticAntigen"
STEP_NAME = {
    "QC": "00-clean",
    "align": "01-star",
    "known": "02-known",
    "novel": "03-novel",
    "salmon": "04-salmon_quant",
    "hla": "05-hla_typing",
    "aeseps": "06-aeSEPs",
    "pvacbind": "07-hla_binding_pred",
}


# ---------------------- Helpers ----------------------
def _script_path(rel_name: str) -> str:
    """
    Resolve the absolute path of a script inside cryptic_pipeline/scripts/.
    Example: _script_path("00-QC.py")
    """
    pkg_path = files("mimicneoai.cryptic_pipeline.scripts")
    return str(pkg_path / rel_name)


def _run_cmd(tool: tools, run_sample_id: str, cmd: List[str], cwd: str | None = None) -> None:
    """
    Execute a shell command via tool.exec_cmd with logging; raise on failure.
    All commands are string-joined safely with shlex.quote.
    """
    cmd = ' '.join(shlex.quote(c) for c in cmd)

    prev_cwd = os.getcwd()
    try:
        if cwd:
            os.chdir(cwd)
        # Single entry point for execution: handles logging, sample ID tagging,
        # stdout/stderr capturing, and returning codes.
        tool.exec_cmd(cmd, run_sample_id, pipline='cryptic')
    except Exception:
        tool.write_log(f"[CMD FAILED]\n{traceback.format_exc()}", "error")
        raise
    finally:
        if cwd:
            os.chdir(prev_cwd)


def _resolve_tumor_control(sample: str, ctrl_from_cfg: str | None) -> Tuple[str, str | None]:
    """
    Resolve tumor and control sample names.
    Supports two formats:
      1) samples list item is "TUMOR,CTRL"
      2) samples list item is "TUMOR", and control comes from configure['others']['control_sample']
    """
    if "," in sample:
        tumor, ctrl = sample.split(",", 1)
        return tumor.strip(), ctrl.strip()
    return sample.strip(), (ctrl_from_cfg.strip() if ctrl_from_cfg else None)


# ---------------------- One-sample pipeline ----------------------
def _run_one_sample(
    sample: str,
    configure: Dict[str, Any],
    paths: Dict[str, Any],
    tool: tools,
) -> None:
    """
    Orchestrate all pipeline steps for a single sample, honoring toggles in the YAML.
    """
    try:
        # ---------- Toggles and core params ----------
        args_cfg = configure["args"]
        others = configure["others"]
        ipt_root = configure["path"]["input_dir"]
        out_root = configure["path"]["output_dir"]

        # Threads (kept simple: using a shared 'threads' entry unless split later)
        n_qc = int(args_cfg.get("threads", 30))
        n_align = int(args_cfg.get("threads", 30))
        n_lncsorf = int(args_cfg.get("threads", 30))
        n_trinity_cpu = int(args_cfg.get("threads", 30))
        n_salmon = int(args_cfg.get("threads", 30))
        n_hla = int(args_cfg.get("threads", 30))
        n_pvacbind = int(args_cfg.get("hla_binding_threads", 5))

        # Feature switches
        do_qc = bool(others.get("QC", True))
        do_align = bool(others.get("alignment", True))
        do_known = bool(others.get("known", True))
        do_novel = bool(others.get("novel", True))
        do_quant = bool(others.get("salmon_quant", True))
        do_quant_ctrl = bool(others.get("salmon_quant_control", True))
        do_hla = bool(others.get("hlatyping", True))
        do_aeseps = bool(others.get("extract_aeseps", True))
        do_pvacbind = bool(others.get("hla_binding_pred", True))

        # Tumor/control resolution
        tumor_sample, ctrl_sample = _resolve_tumor_control(sample, None)

        # External tools (paths)
        TRINITY_SIF = paths["path"]["cryptic"]["TRINITY_SIF"]
        PVACTOOLS_SIF = paths["path"]["common"]["PVACTOOLS"]

        # References (cryptic branch)
        STAR_GENOME_DIR = paths["database"]["cryptic"]["STAR_GENOME_DIR"]

        REF_DIR = paths["database"]["cryptic"]["REF"]["REF_DIR"]
        REF_FA = paths["database"]["cryptic"]["REF"]["REF_FA"]
        REF_GTF = paths["database"]["cryptic"]["REF"]["REF_GTF"]
        REF_LNC_GTF = paths["database"]["cryptic"]["REF"]["REF_LNC_GTF"]
        # Backward compatibility: prefer REF_GENOME if provided, else REF_FA
        REF_GENOME = paths["database"]["cryptic"]["REF"].get("REF_GENOME", REF_FA)

        # HLAHD resources (common)
        HLA_FREQ_DIR = paths["database"]["common"]["HLA"]["FREQ_DATA_DIR"]
        HLA_GENE = paths["database"]["common"]["HLA"]["HLA_GENE"]
        HLA_DICT = paths["database"]["common"]["HLA"]["DICTIONARY"]
        HLA_BOWTIE2_INDEX = paths["database"]["common"]["HLA"]["BOWTIE2_INDEX"]

        # Output layout
        OPT = os.path.join(out_root, tumor_sample)
        DIR00 = os.path.join(OPT, STEP_NAME["QC"])
        DIR01 = os.path.join(OPT, STEP_NAME["align"])
        DIR02_KNOWN = os.path.join(OPT, STEP_NAME["known"])
        DIR03_NOVEL = os.path.join(OPT, STEP_NAME["novel"])
        DIR04 = os.path.join(OPT, STEP_NAME["salmon"])
        DIR05 = os.path.join(OPT, STEP_NAME["hla"])
        DIR06 = os.path.join(OPT, STEP_NAME["aeseps"])
        DIR07 = os.path.join(OPT, STEP_NAME["pvacbind"])
        SHARED = os.path.join(OPT, "023-shared")

        os.makedirs(OPT, exist_ok=True)

        # Raw FASTQs (provided by the caller/configure)
        RAW_R1 = f"{ipt_root}/{tumor_sample}/{tumor_sample}.R1.fq.gz"
        RAW_R2 = f"{ipt_root}/{tumor_sample}/{tumor_sample}.R2.fq.gz"

        RAW_CTRL_R1 = None
        RAW_CTRL_R2 = None
        if ctrl_sample:
            RAW_CTRL_R1 = f"{ipt_root}/{ctrl_sample}/{ctrl_sample}.R1.fq.gz"
            RAW_CTRL_R2 = f"{ipt_root}/{ctrl_sample}/{ctrl_sample}.R2.fq.gz"

        # Derived files
        QC_R1 = os.path.join(DIR00, f"{tumor_sample}.R1.QC.fq.gz")
        QC_R2 = os.path.join(DIR00, f"{tumor_sample}.R2.QC.fq.gz")

        QC_CTRL_R1 = os.path.join(DIR00, f"{ctrl_sample}.R1.QC.fq.gz") if ctrl_sample else None
        QC_CTRL_R2 = os.path.join(DIR00, f"{ctrl_sample}.R2.QC.fq.gz") if ctrl_sample else None

        REF_TX_FA = os.path.join(DIR03_NOVEL, "04a.trace_to_ref", "ref.transcripts.fa")
        CONTIGS_FA = os.path.join(DIR03_NOVEL, "04a.trace_to_ref", "contigs.annot.fa")

        QUANT_TUMOR = os.path.join(DIR04, "salmon_quant", "quant.sf")
        QUANT_CTRL = os.path.join(DIR04, "salmon_quant_control", "quant.sf")
        DEDUP_MAP = os.path.join(DIR04, f"{tumor_sample}.merged_tx.dedup.map.tsv")

        AESEPs_PEP = os.path.join(DIR06, f"{tumor_sample}.aeSEPs.pep")
        ABERRANT_TABLE = os.path.join(DIR06, f"{tumor_sample}.aberrant_noncoding.annot.csv")
        HLA_FINAL_TXT = os.path.join(DIR05, tumor_sample, "result", f"{tumor_sample}_final.result.txt")

        # Minimum requirements for aeSEPs (can be overridden in YAML)
        min_tpm_tumor = float(others.get("min_tpm_tumor", 5.0))
        max_tpm_ctrl = float(others.get("max_tpm_ctrl", 0.5))
        min_log2fc = float(others.get("min_log2fc", 4.0))

        # ---------- 00 QC (tumor) ----------
        if do_qc:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("00-QC.py"),
                "--fq1", RAW_R1,
                "--fq2", RAW_R2,
                "-o", DIR00,
                "-p", str(n_qc),
            ])

        # ---------- 00 QC (control) ----------
        if do_qc and ctrl_sample:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("00-QC.py"),
                "--fq1", RAW_CTRL_R1,
                "--fq2", RAW_CTRL_R2,
                "-o", DIR00,
                "-p", str(n_qc),
            ])

        # ---------- 01 STAR alignment ----------
        if do_align:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("01-alignment.py"),
                "-s", tumor_sample,
                "--clean-dir", DIR00,
                "--genome-dir", STAR_GENOME_DIR,
                "--out-root", DIR01,
                "-p", str(n_align),
            ])

        # Use produced BAM if available (to pass into step 02)
        IN_BAM = os.path.join(DIR01, tumor_sample, f"{tumor_sample}.star", f"{tumor_sample}Aligned.out.bam")
        extra_in_bam = ["--in-bam", IN_BAM] if os.path.isfile(IN_BAM) else []

        # ---------- 02 known lnc/sORF ----------
        if do_known:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("02-lnc_sORF_pipeline.py"),
                "-s", tumor_sample,
                "-m", "known",
                "-o", DIR02_KNOWN,
                "--shared-dir", SHARED,
                "--threads", str(n_lncsorf),
                "--ref-dir", REF_DIR,
                "--ref-fa", REF_FA,
                "--ref-gtf", REF_GTF,
                "--ref-lnc-gtf", REF_LNC_GTF,
                *extra_in_bam,
            ])

        # ---------- 02 novel lnc/sORF (with Trinity) ----------
        if do_novel:
            cmd_novel = [
                sys.executable, _script_path("02-lnc_sORF_pipeline.py"),
                "-s", tumor_sample,
                "-m", "novel",
                "-o", DIR03_NOVEL,
                "--shared-dir", SHARED,
                "--threads", str(n_lncsorf),
                "--trinity-mode", str(others.get("trinity_mode", "apptainer")),
                "--trinity-sif", TRINITY_SIF,
                "--trinity-cpu", str(n_trinity_cpu),
                "--trinity-mem", str(args_cfg.get("trinity_mem", "100G")),
                "--ref-dir", REF_DIR,
                "--ref-fa", REF_FA,
                "--ref-gtf", REF_GTF,
                "--ref-lnc-gtf", REF_LNC_GTF,
                *extra_in_bam,
            ]
            _run_cmd(tool, sample, cmd_novel)

        # ---------- 03 Salmon quant (tumor) ----------
        if do_quant:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("03-salmon_quant.py"),
                "-s", tumor_sample,
                "-o", DIR04,
                "--ref-tx-fa", REF_TX_FA,
                "--contigs-fa", CONTIGS_FA,
                "--genome-fa", REF_GENOME,
                "--fq1", QC_R1,
                "--fq2", QC_R2,
                "--threads", str(n_salmon),
                "--kmer", str(int(others.get("salmon_kmer", 31))),
            ])

        # ---------- 04 Salmon quant (control) ----------
        if do_quant_ctrl and ctrl_sample:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("04-salmon_quant_control.py"),
                "--fq1", QC_CTRL_R1,
                "--fq2", QC_CTRL_R2,
                "-i", os.path.join(DIR04, "salmon_index"),
                "-o", os.path.join(DIR04, "salmon_quant_control"),
                "-p", str(n_salmon),
            ])

        # ---------- 05 HLA typing (hlahd) ----------
        if do_hla:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("05-hla_typing.py"),
                "-s", tumor_sample,
                "--r1", QC_R1,
                "--r2", QC_R2,
                "--output-dir", OPT,
                "--step-name-hla", os.path.basename(DIR05),
                "-t", str(n_hla),
                "--freq-data-dir", HLA_FREQ_DIR,
                "--HLA-gene", HLA_GENE,
                "--dictionary", HLA_DICT,
                "--hla-gen", HLA_BOWTIE2_INDEX,
            ])

        # ---------- 06 Extract aeSEPs ----------
        if do_aeseps:
            _run_cmd(tool, sample, [
                sys.executable, _script_path("06-sORF-encoded_peptides.py"),
                "--quant-tumor", QUANT_TUMOR,
                "--quant-ctrl", QUANT_CTRL if ctrl_sample else "",
                "--dedup-map", DEDUP_MAP,
                "--known-tx-fa", os.path.join(DIR02_KNOWN, "04k.bcf_consensus", "lncRNA.consensus.fa"),
                "--pep-known", os.path.join(DIR02_KNOWN, "06.filter_sORFs", f"{tumor_sample}.lncRNA.sORFs.pep"),
                "--pep-novel", os.path.join(DIR03_NOVEL, "06.filter_sORFs", f"{tumor_sample}.lncRNA.sORFs.pep"),
                "--out-dir", DIR06,
                "--out-name", os.path.basename(AESEPs_PEP),
                "--save-table", os.path.basename(ABERRANT_TABLE),
                "--min-tpm-tumor", str(min_tpm_tumor),
                "--max-tpm-ctrl", str(max_tpm_ctrl),
                "--min-log2fc", str(min_log2fc),
            ])

        # ---------- 07 HLA binding prediction (pvacbind) ----------
        if do_pvacbind:
            algos = others.get(
                "algo",
                "BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI "
                "MHCnuggetsII NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan "
                "NetMHCpanEL PickPocket SMM SMMPMBEC",
            )
            e1_lengths = others.get("mhcI_lengths", "8,9,10")
            e2_lengths = others.get("mhcII_lengths", "15")

            _run_cmd(tool, sample, [
                sys.executable, _script_path("07-hla_binding_pred.py"),
                "-s", tumor_sample,
                "--pep-fasta", AESEPs_PEP,
                "--hla-file", HLA_FINAL_TXT,
                "-o", DIR07,
                "--pvactools", PVACTOOLS_SIF,
                "-t", str(n_pvacbind),
                "--algos", algos,
                "--e1-lengths", e1_lengths,
                "--e2-lengths", e2_lengths,
            ])

        tool.write_log(f"[DONE] Completed cryptic pipeline: {OPT}", "info")

    except Exception:
        tool.write_log(f"Worker error for sample '{sample}':\n{traceback.format_exc()}", "error")
        raise


# -------- Non-daemon Pool (optional: needed only if you spawn pools inside workers) --------
class _NoDaemonProcess(multiprocessing.Process):
    """A Process whose daemon attribute is always False (allows child processes)."""
    # For Python < 3.8
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


if hasattr(multiprocessing, "get_start_method"):
    # Python ≥ 3.8: override Process constructor on Pool to use _NoDaemonProcess
    class NoDaemonPool(multiprocessing.pool.Pool):
        @staticmethod
        def Process(_, *args, **kwargs):
            return _NoDaemonProcess(*args, **kwargs)
else:
    # Python < 3.8 fallback
    class NoDaemonPool(multiprocessing.pool.Pool):
        Process = _NoDaemonProcess


# ---------------------- Multiprocess runner ----------------------
def _run_pool(samples: List[str], pool_size: int, configure, paths, tool: tools) -> None:
    """
    Submit one asynchronous task per sample and wait for completion.
    Errors inside workers are reported via tool.print_pool_error.
    """
    with NoDaemonPool(processes=pool_size) as pool:
        for s in samples:
            pool.apply_async(
                _run_one_sample,
                (s, configure, paths, tool),
                error_callback=tool.print_pool_error,
            )
        pool.close()
        pool.join()


def _peek_output_dir(cfg_path: str) -> str | None:
    """
    Quickly peek 'configure' YAML to obtain path.output_dir for choosing a working directory.
    Return None on any failure; do not raise.
    """
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        out_dir = (data.get("path") or {}).get("output_dir")
        if not out_dir:
            return None
        return str(Path(out_dir).expanduser())
    except Exception:
        return None


# ---------------------- CLI entry ----------------------
def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="MimicNeoAI (Cryptic pipeline)")
    parser.add_argument("-c", "--configure", required=True, help="Path to configuration YAML")
    parser.add_argument("-p", "--paths", dest="paths", required=False, help="Path to paths YAML")
    parser.add_argument("--workdir", default=None, help="Working directory (default: config.path.output_dir or CWD)")
    args = parser.parse_args(argv)

    # 1) Decide working/log root: --workdir > configure.path.output_dir > CWD
    cfg_output_dir = _peek_output_dir(args.configure)
    workdir = args.workdir or cfg_output_dir or os.getcwd()
    workdir = str(Path(workdir).resolve())

    # Manager + tool/logger (log root = workdir)
    mgr = Manager()
    lock = mgr.Lock()
    tool_obj = tools(workdir, FLAG, lock)

    tool_obj.write_log(f"work_path: {tool_obj.sys_path}", "info")
    tool_obj.write_log(f"start_log: {tool_obj.start_log}", "info")
    tool_obj.write_log(f"cmd_log: {tool_obj.cmd_log_dir}", "info")
    if args.workdir is None and cfg_output_dir:
        tool_obj.write_log(f"workdir inferred from configure.path.output_dir: {cfg_output_dir}", "info")

    # Load YAMLs
    configure = tool_obj.get_configure(args.configure)
    paths = tool_obj.get_paths(args.paths)
    tool_obj.write_log(f"configures: {configure}", "info")
    tool_obj.write_log(f"paths: {paths}", "info")

    # Samples & pool
    samples = [str(s) for s in configure["samples"]]
    pool_size = int(configure.get("args", {}).get("pool_size", 1))

    # Inject step names (non-destructive)
    configure = dict(configure)
    configure.setdefault("step_name", {}).update(STEP_NAME)

    # Share variables across processes
    tool_obj.sharing_variable(mgr, samples)

    # Run
    _run_pool(samples, pool_size, configure, paths, tool_obj)

    tool_obj.summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
