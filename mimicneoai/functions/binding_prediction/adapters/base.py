"""Base classes and common helpers for binding predictor adapters."""

from __future__ import annotations

import csv
import gzip
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

from mimicneoai.functions.binding_prediction.schema import (
    BindingPrediction,
    PredictionJob,
    write_prediction_rows,
)


@dataclass(frozen=True)
class AdapterConfig:
    """Filesystem configuration for local predictor tools."""

    mhcflurry_predict_bin: str = "/workspace/pkgs/mhcflurry/.venv/bin/mhcflurry-predict"
    mhcnuggets_python_bin: str = "/workspace/pkgs/mhcnuggets/.venv/bin/python"
    mhcnuggets_script: str = "/workspace/pkgs/mhcnuggets/mhcnuggets/src/predict.py"
    mhcnuggets_cwd: str = "/workspace/pkgs/mhcnuggets"
    mhcnuggets_rank_output: bool = True
    netmhcpan_bin: str = "/workspace/pkgs/NetMHCpan/netMHCpan-4.2/netMHCpan"
    netmhciipan_bin: str = "/workspace/pkgs/NetMHCIIpan/netMHCIIpan-4.3/netMHCIIpan"
    iedb_mhci_script: str = "/workspace/pkgs/IEDB/mhc_i/src/predict_binding.py"
    iedb_mhci_cwd: str = "/workspace/pkgs/IEDB/mhc_i"
    iedb_mhcii_script: str = "/workspace/pkgs/IEDB/mhc_ii/mhc_II_binding.py"
    iedb_mhcii_cwd: str = "/workspace/pkgs/IEDB/mhc_ii"
    python_bin: str = "python"
    iedb_mhcii_python_bin: str = "/workspace/pkgs/mhcflurry/.venv/bin/python"
    chunk_size: Optional[int] = None
    resume: bool = True
    device: str = "cpu"
    gpu_id: str = "0"
    tf_intra_op_threads: int = 1
    tf_inter_op_threads: int = 1
    command_timeout: Optional[int] = None


class PredictorAdapter:
    """Base interface for one family of prediction algorithms."""

    supported_algorithms: frozenset[str] = frozenset()
    default_chunk_size: int = 5000

    def __init__(self, config: AdapterConfig):
        self.config = config

    def run_job(self, job: PredictionJob) -> Path:
        raise NotImplementedError

    def chunk_size(self) -> int:
        return self.config.chunk_size or self.default_chunk_size


def predictor_env(config: Optional[AdapterConfig] = None, extra: Optional[dict[str, str]] = None) -> dict[str, str]:
    """Return an environment for predictor subprocesses.

    CPU is the default. GPU use must be explicitly requested through
    ``AdapterConfig.device == "gpu"`` so clinical batch runs remain portable on
    servers without CUDA devices.
    """

    env = os.environ.copy()
    device = (config.device if config else "cpu").lower()
    if device == "gpu":
        cuda_visible_devices = config.gpu_id if config else "0"
    elif device == "cpu":
        cuda_visible_devices = ""
    else:
        raise ValueError(f"Unsupported predictor device: {device}. Expected cpu or gpu.")
    env.update(
        {
            "CUDA_VISIBLE_DEVICES": cuda_visible_devices,
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "TF_NUM_INTRAOP_THREADS": str(config.tf_intra_op_threads if config else 1),
            "TF_NUM_INTEROP_THREADS": str(config.tf_inter_op_threads if config else 1),
        }
    )
    if extra:
        env.update(extra)
    return env


def ensure_job_dirs(job: PredictionJob) -> None:
    for path in (job.raw_path, job.normalized_path, job.log_path, job.input_path):
        path.parent.mkdir(parents=True, exist_ok=True)


def write_peptide_lines(path: Path, peptides: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for peptide in peptides:
            handle.write(f"{peptide}\n")


def write_peptide_fasta(path: Path, peptides: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for index, peptide in enumerate(peptides, start=1):
            handle.write(f">pep_{index:06d}\n{peptide}\n")


def run_logged_command(
    cmd: Sequence[str],
    cwd: Optional[Path],
    log_path: Path,
    stdout_path: Optional[Path] = None,
    config: Optional[AdapterConfig] = None,
) -> None:
    """Run a subprocess, capturing command, stdout, and stderr."""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        log_handle.write("CMD: " + " ".join(cmd) + "\n")
        log_handle.flush()
        if stdout_path is None:
            proc = subprocess.run(
                list(cmd),
                cwd=str(cwd) if cwd else None,
                env=predictor_env(config),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=config.command_timeout if config else None,
            )
        else:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with stdout_path.open("w") as stdout_handle:
                proc = subprocess.run(
                    list(cmd),
                    cwd=str(cwd) if cwd else None,
                    env=predictor_env(config),
                    stdout=stdout_handle,
                    stderr=log_handle,
                    text=True,
                    check=False,
                    timeout=config.command_timeout if config else None,
                )
        if proc.returncode != 0:
            raise RuntimeError(f"command failed with exit code {proc.returncode}: {' '.join(cmd)}")


def read_tabular_after_header(path: Path, header_prefix: str) -> list[dict[str, str]]:
    """Read a tabular output file after the first header line matching a prefix."""

    lines = path.read_text().splitlines()
    header_index = None
    for index, line in enumerate(lines):
        if line.startswith(header_prefix):
            header_index = index
            break
    if header_index is None:
        return []
    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    return list(reader)


def parse_netmhc_stdout_table(path: Path) -> list[dict[str, str]]:
    """Parse whitespace-aligned NetMHC stdout tables.

    NetMHCpan and NetMHCIIpan write human-readable tables to stdout. The final
    BindLevel field can contain spaces, for example "<= WB"; this parser follows
    the same rule as the historical merge scripts: map the first n-1 tokens to
    the first n-1 header fields, then join all remaining tokens into the last
    field.
    """

    header: Optional[list[str]] = None
    rows: list[dict[str, str]] = []
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("Pos"):
                header = stripped.split()
                continue
            if header is None:
                continue
            if not stripped or stripped.startswith("-") or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if not parts:
                continue
            try:
                int(parts[0])
            except ValueError:
                continue
            column_count = len(header)
            if len(parts) < column_count:
                parts = parts + [""] * (column_count - len(parts))
                row = parts[:column_count]
            else:
                row = parts[: column_count - 1] + [" ".join(parts[column_count - 1 :])]
            rows.append(dict(zip(header, row)))
    return rows


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def write_error_predictions(job: PredictionJob, error: Exception) -> Path:
    predictions = [
        BindingPrediction(
            peptide=peptide,
            hla_allele=job.hla_allele,
            algorithm=algorithm,
            mhc_class=job.mhc_class,
            peptide_length=job.peptide_length,
            status="error",
            raw_file=str(job.raw_path),
            error=str(error),
        )
        for peptide in job.peptides
        for algorithm in (job.output_algorithms or (job.algorithm,))
    ]
    write_prediction_rows(job.normalized_path, predictions)
    return job.normalized_path


def write_ok_predictions(job: PredictionJob, predictions: Iterable[BindingPrediction]) -> Path:
    write_prediction_rows(job.normalized_path, predictions)
    return job.normalized_path
