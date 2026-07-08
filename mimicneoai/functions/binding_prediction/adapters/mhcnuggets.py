"""MHCnuggets adapter for class I and class II affinity prediction."""

from __future__ import annotations

import csv
from pathlib import Path

from mimicneoai.functions.binding_prediction.schema import BindingPrediction, PredictionJob

from .base import (
    PredictorAdapter,
    ensure_job_dirs,
    run_logged_command,
    write_error_predictions,
    write_ok_predictions,
    write_peptide_lines,
)


CLASS_BY_ALGORITHM = {
    "MHCnuggetsI": "I",
    "MHCnuggetsII": "II",
}


class MhcnuggetsAdapter(PredictorAdapter):
    """Run MHCnuggets using binding-affinity models."""

    supported_algorithms = frozenset(CLASS_BY_ALGORITHM)
    default_chunk_size = 3000
    class_i_min_length = 8
    class_i_max_length = 15

    def run_job(self, job: PredictionJob) -> Path:
        if self.config.resume and job.normalized_path.exists() and job.normalized_path.stat().st_size > 0:
            return job.normalized_path
        ensure_job_dirs(job)
        try:
            expected_class = "MHC-I" if job.algorithm == "MHCnuggetsI" else "MHC-II"
            if job.mhc_class != expected_class:
                raise ValueError(f"{job.algorithm} expects {expected_class}, got {job.mhc_class}")
            if job.algorithm == "MHCnuggetsI" and not (
                self.class_i_min_length <= job.peptide_length <= self.class_i_max_length
            ):
                raise ValueError(
                    f"MHCnuggetsI supports class-I peptide length "
                    f"{self.class_i_min_length}-{self.class_i_max_length}, got {job.peptide_length}"
                )
            write_peptide_lines(job.input_path, job.peptides)
            rank_error = ""
            if self.config.mhcnuggets_rank_output:
                try:
                    run_logged_command(
                        build_mhcnuggets_command(job, rank_output=True, config=self.config),
                        Path(self.config.mhcnuggets_cwd),
                        job.log_path,
                        config=self.config,
                    )
                except Exception as exc:  # noqa: BLE001
                    rank_error = str(exc)
                    run_logged_command(
                        build_mhcnuggets_command(job, rank_output=False, config=self.config),
                        Path(self.config.mhcnuggets_cwd),
                        job.log_path,
                        config=self.config,
                    )
            else:
                run_logged_command(
                    build_mhcnuggets_command(job, rank_output=False, config=self.config),
                    Path(self.config.mhcnuggets_cwd),
                    job.log_path,
                    config=self.config,
                )
            rows = list(csv.DictReader(job.raw_path.open()))
            rank_by_peptide = read_mhcnuggets_rank_output(job.raw_path) if self.config.mhcnuggets_rank_output else {}
            rows_by_peptide: dict[str, dict[str, str]] = {}
            for row in rows:
                peptide = normalize_peptide(row.get("peptide", ""))
                if peptide and peptide not in rows_by_peptide:
                    rows_by_peptide[peptide] = row

            predictions = []
            missing_peptides = []
            for peptide in job.peptides:
                normalized_peptide = normalize_peptide(peptide)
                row = rows_by_peptide.get(normalized_peptide)
                if row is None:
                    missing_peptides.append(peptide)
                    continue
                rank = rank_by_peptide.get(normalized_peptide, "")
                predictions.append(
                    BindingPrediction(
                        peptide=peptide,
                        hla_allele=job.hla_allele,
                        algorithm=job.algorithm,
                        mhc_class=job.mhc_class,
                        peptide_length=job.peptide_length,
                        ic50=row.get("ic50", ""),
                        ic50_type="MHCnuggets BA IC50(nM)",
                        percentile=rank,
                        percentile_type="MHCnuggets human_proteome_rank" if rank else "",
                        raw_rank=rank,
                        raw_rank_type="MHCnuggets human_proteome_rank" if rank else "",
                        status="partial_ok" if rank_error else "ok",
                        raw_file=str(job.raw_path),
                        error=f"MHCnuggets rank output failed; IC50-only fallback used: {rank_error}"
                        if rank_error
                        else "",
                    )
                )
            for peptide in missing_peptides:
                predictions.append(
                    BindingPrediction(
                        peptide=peptide,
                        hla_allele=job.hla_allele,
                        algorithm=job.algorithm,
                        mhc_class=job.mhc_class,
                        peptide_length=job.peptide_length,
                        status="error",
                        raw_file=str(job.raw_path),
                        error="MHCnuggets did not return score for this peptide",
                    )
                )
            return write_ok_predictions(job, predictions)
        except Exception as exc:  # noqa: BLE001
            return write_error_predictions(job, exc)


def normalize_peptide(peptide: str) -> str:
    return str(peptide or "").strip().upper()


def build_mhcnuggets_command(job: PredictionJob, rank_output: bool, config) -> list[str]:
    cmd = [
        config.mhcnuggets_python_bin,
        config.mhcnuggets_script,
        "-c",
        CLASS_BY_ALGORITHM[job.algorithm],
        "-p",
        str(job.input_path),
        "-a",
        format_mhcnuggets_allele(job.hla_allele),
        "-o",
        str(job.raw_path),
        "-M",
    ]
    if rank_output:
        cmd.extend(["-r", "True"])
    return cmd


def format_mhcnuggets_allele(allele: str) -> str:
    """Convert common HLA notation to MHCnuggets model naming."""

    parts = []
    for part in allele.strip().split("/"):
        part = part.strip()
        if part.startswith("HLA-"):
            part = part[4:]
        parts.append("HLA-" + part.replace("*", ""))
    return "-".join(parts)


def read_mhcnuggets_rank_output(output_path: Path) -> dict[str, str]:
    """Read the optional MHCnuggets rank file matching an output path."""

    rank_path = first_existing_mhcnuggets_rank_path(output_path)
    if not rank_path.exists():
        return {}
    result: dict[str, str] = {}
    with rank_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            peptide = normalize_peptide(row.get("peptide", ""))
            if peptide:
                result[peptide] = row.get("human_proteome_rank", "")
    return result


def mhcnuggets_rank_path(output_path: Path) -> Path:
    if output_path.suffix:
        return output_path.with_name(f"{output_path.stem}_ranks{output_path.suffix}")
    return Path(str(output_path) + "_ranks")


def mhcnuggets_actual_rank_path(output_path: Path) -> Path:
    """Return the rank filename produced by upstream MHCnuggets.

    The upstream script uses "{}_ranks.{}".format(file_name, extension), where
    extension already includes a leading dot. For an output named
    chunk.raw.tsv this yields chunk.raw_ranks..tsv.
    """

    if output_path.suffix:
        return output_path.with_name(f"{output_path.stem}_ranks.{output_path.suffix}")
    return Path(str(output_path) + "_ranks")


def first_existing_mhcnuggets_rank_path(output_path: Path) -> Path:
    for candidate in (mhcnuggets_rank_path(output_path), mhcnuggets_actual_rank_path(output_path)):
        if candidate.exists():
            return candidate
    return mhcnuggets_rank_path(output_path)
