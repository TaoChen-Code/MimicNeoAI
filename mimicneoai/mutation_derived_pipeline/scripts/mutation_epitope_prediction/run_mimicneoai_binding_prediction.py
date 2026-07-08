"""Run the MimicNeoAI mutation-derived binding prediction workflow.

The workflow uses pVACtools only as an external source generator for mutation
annotation and WT/MT protein FASTA. MimicNeoAI then builds mutation-covering
epitope windows, runs local binding predictors, and writes a pVACseq-compatible
merged all_epitopes table.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional


DEFAULT_ALGORITHMS = (
    "MHCflurry",
    "MHCflurryEL",
    "MHCnuggetsI",
    "MHCnuggetsII",
    "NNalign",
    "NetMHCpan",
    "NetMHCpanEL",
    "NetMHCIIpan",
    "NetMHCIIpanEL",
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[4]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", required=True)
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--hla-file", required=True)
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Workflow output directory, e.g. <sample>/07.binding_prediction_mimicneoai.",
    )
    parser.add_argument("--pvactools-sif", required=True)
    parser.add_argument("--apptainer", default="apptainer")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--bind", action="append", default=[])
    parser.add_argument("--mhc-i-lengths", default="8,9,10,11")
    parser.add_argument("--mhc-ii-lengths", default="15")
    parser.add_argument("--protein-flank-length", type=int, default=25)
    parser.add_argument("--extended-length", type=int, default=27)
    parser.add_argument("--algorithms", default=",".join(DEFAULT_ALGORITHMS))
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--mt-workers", type=int, default=None)
    parser.add_argument("--wt-workers", type=int, default=None)
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--command-timeout", type=int, default=None)
    parser.add_argument("--netmhcpan-bin", default="netMHCpan")
    parser.add_argument("--netmhciipan-bin", default="netMHCIIpan")
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--pvacseq-merged", default=None, help="Optional legacy pVACseq merged TSV for validation.")
    parser.add_argument("--no-pass-only", action="store_true")
    parser.add_argument("--no-archive-runner-workdirs", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = Path(args.outdir)
    started_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    paths = workflow_paths(outdir, args.sample, args.protein_flank_length)
    commands = build_commands(args, paths)
    if args.dry_run:
        for label, command in commands:
            print(f"[DRY-RUN] {label}: {format_command(command)}", flush=True)
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    for key in ("epitope_tasks", "binding_predictions", "merged_epitopes", "archive"):
        paths[key].mkdir(parents=True, exist_ok=True)

    run_steps(commands, paths, python_bin=args.python_bin)
    archive_runner_workdirs(paths, args.sample, enabled=not args.no_archive_runner_workdirs)
    summary = build_workflow_summary(args, paths, started_at)
    summary_path = outdir / f"{args.sample}.binding_prediction_mimicneoai.summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(f"[DONE] MimicNeoAI binding prediction workflow finished: {outdir}", flush=True)
    print(f"[DONE] Summary: {summary_path}", flush=True)
    return 0


def workflow_paths(outdir: Path, sample: str, protein_flank_length: int) -> dict[str, Path]:
    return {
        "root": outdir,
        "pvactools_sources": outdir / "01_pvactools_sources",
        "epitope_tasks": outdir / "02_epitope_tasks",
        "binding_predictions": outdir / "03_binding_predictions",
        "merged_epitopes": outdir / "04_merged_epitopes",
        "archive": outdir / "archive",
        "converter_tsv": outdir / "01_pvactools_sources" / f"{sample}.pvacseq_converter.tsv",
        "protein_fasta": outdir / "01_pvactools_sources" / f"{sample}.protein.flank{protein_flank_length}.wt_mt.fasta",
        "binding_tasks": outdir / "02_epitope_tasks" / "binding_tasks.tsv",
        "epitope_windows": outdir / "02_epitope_tasks" / "epitope_windows.tsv",
        "variant_annotation": outdir / "02_epitope_tasks" / "variant_annotation.tsv",
        "mt_tasks": outdir / "03_binding_predictions" / f"{sample}.MT.binding_tasks.tsv",
        "wt_tasks": outdir / "03_binding_predictions" / f"{sample}.WT.binding_tasks.tsv",
        "mt_runner": outdir / "03_binding_predictions" / f"{sample}.MT.runner",
        "wt_runner": outdir / "03_binding_predictions" / f"{sample}.WT.runner",
        "mt_predictions": outdir / "03_binding_predictions" / f"{sample}.MT.binding_predictions.long.tsv",
        "wt_predictions": outdir / "03_binding_predictions" / f"{sample}.WT.binding_predictions.long.tsv",
        "merged": outdir / "04_merged_epitopes" / f"{sample}.merged.all_epitopes.tsv",
    }


def build_commands(args: argparse.Namespace, paths: dict[str, Path]) -> list[tuple[str, list[str]]]:
    python_bin = args.python_bin
    bind_args: list[str] = []
    for bind_path in args.bind:
        bind_args.extend(["--bind", bind_path])
    pass_args = ["--no-pass-only"] if args.no_pass_only else []
    chunk_args = ["--chunk-size", str(args.chunk_size)] if args.chunk_size else []
    timeout_args = ["--command-timeout", str(args.command_timeout)] if args.command_timeout else []
    pvacseq_validation_args = ["--pvacseq-merged", args.pvacseq_merged] if args.pvacseq_merged else []

    commands = [
        (
            "00_prepare_pvacseq_sources",
            [
                python_bin,
                str(SCRIPT_DIR / "00_prepare_pvacseq_sources.py"),
                "-s",
                args.sample,
                "--input-vcf",
                args.input_vcf,
                "-o",
                str(paths["root"]),
                "--pvactools-sif",
                args.pvactools_sif,
                "--apptainer",
                args.apptainer,
                "--bcftools",
                args.bcftools,
                "--tabix",
                args.tabix,
                "--flank-length",
                str(args.protein_flank_length),
                *bind_args,
                *pass_args,
            ],
        ),
        (
            "01_build_epitope_tasks",
            [
                python_bin,
                str(SCRIPT_DIR / "01_build_epitope_tasks.py"),
                "-s",
                args.sample,
                "--converter-tsv",
                str(paths["converter_tsv"]),
                "--protein-fasta",
                str(paths["protein_fasta"]),
                "--hla-file",
                args.hla_file,
                "-o",
                str(paths["epitope_tasks"]),
                "--mhc-i-lengths",
                args.mhc_i_lengths,
                "--mhc-ii-lengths",
                args.mhc_ii_lengths,
                "--protein-flank-length",
                str(args.protein_flank_length),
                "--extended-length",
                str(args.extended_length),
                "--algorithms",
                args.algorithms,
                *pvacseq_validation_args,
            ],
        ),
        (
            "02_split_binding_tasks",
            [
                python_bin,
                str(SCRIPT_DIR / "02_split_binding_tasks.py"),
                "-s",
                args.sample,
                "--binding-tasks",
                str(paths["binding_tasks"]),
                "--epitope-windows",
                str(paths["epitope_windows"]),
                "-o",
                str(paths["binding_predictions"]),
                "--mt-output",
                str(paths["mt_tasks"]),
                "--wt-output",
                str(paths["wt_tasks"]),
            ],
        ),
        (
            "03_predict_MT",
            [
                python_bin,
                "-m",
                "mimicneoai.functions.binding_prediction.runner",
                "--tasks",
                str(paths["mt_tasks"]),
                "-o",
                str(paths["mt_runner"]),
                "--workers",
                str(args.mt_workers or args.workers),
                "--device",
                args.device,
                "--netmhcpan-bin",
                args.netmhcpan_bin,
                "--netmhciipan-bin",
                args.netmhciipan_bin,
                *chunk_args,
                *timeout_args,
            ],
        ),
        (
            "03_predict_WT",
            [
                python_bin,
                "-m",
                "mimicneoai.functions.binding_prediction.runner",
                "--tasks",
                str(paths["wt_tasks"]),
                "-o",
                str(paths["wt_runner"]),
                "--workers",
                str(args.wt_workers or args.workers),
                "--device",
                args.device,
                "--netmhcpan-bin",
                args.netmhcpan_bin,
                "--netmhciipan-bin",
                args.netmhciipan_bin,
                *chunk_args,
                *timeout_args,
            ],
        ),
        (
            "04_merge_binding_predictions",
            [
                python_bin,
                str(SCRIPT_DIR / "02_merge_binding_predictions.py"),
                "--variant-annotation",
                str(paths["variant_annotation"]),
                "--epitope-windows",
                str(paths["epitope_windows"]),
                "--binding-predictions",
                str(paths["mt_predictions"]),
                "--binding-predictions",
                str(paths["wt_predictions"]),
                "--binding-tasks",
                str(paths["mt_tasks"]),
                "--output-profile",
                "pvacseq",
                "-o",
                str(paths["merged"]),
                "--summary",
                str(paths["merged_epitopes"] / f"{args.sample}.merged.summary.json"),
                *pvacseq_validation_args,
            ],
        ),
    ]
    return commands


def run_steps(commands: list[tuple[str, list[str]]], paths: dict[str, Path], python_bin: str) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    for label, command in commands:
        print(f"[RUN] {label}: {format_command(command)}", flush=True)
        subprocess.run(command, check=True, env=env)
        if label == "03_predict_MT":
            move_runner_outputs(paths["mt_runner"], paths["binding_predictions"], "MT")
        elif label == "03_predict_WT":
            move_runner_outputs(paths["wt_runner"], paths["binding_predictions"], "WT")


def move_runner_outputs(runner_dir: Path, parent: Path, source_type: str) -> None:
    sample = runner_dir.name.rsplit(f".{source_type}.runner", 1)[0]
    mapping = {
        runner_dir / "binding_predictions.long.tsv": parent / f"{sample}.{source_type}.binding_predictions.long.tsv",
        runner_dir / "binding_predictions.summary.json": parent / f"{sample}.{source_type}.binding_predictions.summary.json",
    }
    for source, target in mapping.items():
        if not source.exists():
            raise FileNotFoundError(f"Expected runner output missing: {source}")
        if target.exists():
            target.unlink()
        source.replace(target)


def archive_runner_workdirs(paths: dict[str, Path], sample: str, enabled: bool) -> None:
    if not enabled:
        return
    for source_type in ("MT", "WT"):
        runner_dir = paths[f"{source_type.lower()}_runner"]
        if not runner_dir.exists():
            continue
        target = paths["archive"] / "03_binding_predictions" / runner_dir.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            target = target.with_name(f"{target.name}.{timestamp}.old")
        shutil.move(str(runner_dir), str(target))


def build_workflow_summary(args: argparse.Namespace, paths: dict[str, Path], started_at: str) -> dict[str, object]:
    mt_split_summary = read_json(paths["binding_predictions"] / f"{args.sample}.split_binding_tasks.summary.json")
    mt_prediction_summary = read_json(paths["binding_predictions"] / f"{args.sample}.MT.binding_predictions.summary.json")
    wt_prediction_summary = read_json(paths["binding_predictions"] / f"{args.sample}.WT.binding_predictions.summary.json")
    merged_summary = read_json(paths["merged_epitopes"] / f"{args.sample}.merged.summary.json")
    return {
        "sample": args.sample,
        "started_at": started_at,
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "workflow": "MimicNeoAI-native mutation-derived binding prediction",
        "output_dir": str(paths["root"]),
        "pvactools_boundary": "pVACtools is used only for VCF annotation conversion and WT/MT protein FASTA source generation.",
        "pvactools_sif": args.pvactools_sif,
        "pvactools_version": capture_command([args.apptainer, "exec", str(args.pvactools_sif), "pvacseq", "--version"]),
        "netmhcpan_version": capture_command([args.netmhcpan_bin]),
        "netmhciipan_version": capture_command([args.netmhciipan_bin]),
        "algorithms": [item for item in args.algorithms.replace(",", " ").split() if item],
        "mhc_i_lengths": args.mhc_i_lengths,
        "mhc_ii_lengths": args.mhc_ii_lengths,
        "workers": args.workers,
        "mt_workers": args.mt_workers or args.workers,
        "wt_workers": args.wt_workers or args.workers,
        "mt_task_rows": mt_split_summary.get("mt_task_rows"),
        "wt_task_rows": mt_split_summary.get("wt_task_rows"),
        "mt_prediction_rows": mt_prediction_summary.get("prediction_rows"),
        "wt_prediction_rows": wt_prediction_summary.get("prediction_rows"),
        "merged_rows": merged_summary.get("rows"),
        "fs_wt_rule": "Frameshift WT Epitope Seq and corresponding WT IC50/fold-change/percentile fields are left blank.",
        "outputs": {
            "mt_binding_predictions": str(paths["mt_predictions"]),
            "wt_binding_predictions": str(paths["wt_predictions"]),
            "merged_all_epitopes": str(paths["merged"]),
        },
    }


def read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    with path.open() as handle:
        return json.load(handle)


def capture_command(command: list[str]) -> str:
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=30, check=False)
    except Exception as exc:  # pragma: no cover - best-effort metadata only
        return f"unavailable: {exc}"
    text = "\n".join(part.strip() for part in (result.stdout, result.stderr) if part.strip())
    return text[:2000]


def format_command(command: list[str]) -> str:
    return " ".join(str(part) for part in command)


if __name__ == "__main__":
    raise SystemExit(main())
