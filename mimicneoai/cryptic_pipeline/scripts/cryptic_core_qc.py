#!/usr/bin/env python3
"""Materialize a strict, pre-binding Cryptic Core.

This stage consumes aeSEP parent ORFs after ORF-genome filtering and writes
pre-tiled HLA-I/HLA-II peptide Core FASTA files. It intentionally keeps RNA
editing, junction support, and external atlas checks as explicit not-evaluated
fields until those resources are configured.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


POLICY_VERSION = "cryptic_core_qc_v1.0"
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
ACCEPTED_SOURCE_CLASSES = {"novel", "noncoding"}


def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def parse_int_list(value: str) -> list[int]:
    out: list[int] = []
    for token in str(value).replace(",", " ").split():
        token = token.strip()
        if token:
            out.append(int(token))
    if not out:
        raise ValueError("At least one peptide length is required")
    return sorted(set(out))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"path": str(path), "exists": False}
    return {
        "path": str(path),
        "exists": True,
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def read_fasta(path: Path) -> Iterable[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    header: Optional[str] = None
    chunks: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def fasta_record_id(header: str) -> str:
    fields = str(header).strip().split()
    if not fields:
        raise ValueError("FASTA header does not contain a record id")
    return fields[0]


def source_parent_id(record_id: str) -> str:
    return record_id.rsplit(".p", 1)[0] if ".p" in record_id else record_id


def normalize_parent_sequence(sequence: str) -> tuple[str, list[str]]:
    seq = str(sequence).replace(" ", "").replace("\n", "").strip().upper()
    reasons: list[str] = []
    if not seq:
        return "", ["empty_sequence"]
    if seq.endswith("*"):
        seq = seq[:-1]
    if "*" in seq:
        reasons.append("internal_stop")
    if "-" in seq:
        reasons.append("gap_character")
    bad = sorted(set(seq).difference(CANONICAL_AA))
    if bad:
        reasons.append("noncanonical_amino_acid:" + ",".join(bad))
    if not seq:
        reasons.append("empty_after_terminal_stop_trim")
    return seq, reasons


def iter_windows(sequence: str, lengths: Iterable[int]) -> Iterable[tuple[int, int, str]]:
    for length in lengths:
        if len(sequence) < length:
            continue
        for start in range(0, len(sequence) - length + 1):
            yield start + 1, length, sequence[start : start + length]


def load_human_reference_matches(
    human_proteome_fasta: Optional[Path],
    candidate_peptides: set[str],
) -> tuple[set[str], dict[str, object]]:
    if not human_proteome_fasta:
        return set(), {"status": "not_evaluated", "path": ""}
    if not human_proteome_fasta.exists():
        raise FileNotFoundError(f"human reference proteome FASTA not found: {human_proteome_fasta}")
    matches: set[str] = set()
    lengths = sorted({len(peptide) for peptide in candidate_peptides})
    records = 0
    for _header, raw_seq in read_fasta(human_proteome_fasta):
        records += 1
        seq, reasons = normalize_parent_sequence(raw_seq)
        if reasons or not seq:
            continue
        for length in lengths:
            if len(seq) < length:
                continue
            for start in range(0, len(seq) - length + 1):
                pep = seq[start : start + length]
                if pep in candidate_peptides:
                    matches.add(pep)
    return matches, {
        "status": "evaluated",
        "path": str(human_proteome_fasta),
        "sha256": sha256_file(human_proteome_fasta),
        "records_scanned": records,
        "candidate_peptides": len(candidate_peptides),
        "matched_peptides": len(matches),
    }


def write_fasta(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(f">{row['parent_record_id']}\n{row['peptide']}\n")


def write_peptide_core_fasta(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(f">{row['peptide_record_id']}\n{row['peptide']}\n")


def write_tsv(rows: list[dict[str, object]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_json(path: Path) -> dict[str, object]:
    with path.open() as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def output_signature(paths: list[Path]) -> dict[str, dict[str, object]]:
    return {path.name: file_identity(path) for path in paths}


def outputs_match_signature(signature: dict[str, object]) -> bool:
    for identity in signature.values():
        if not isinstance(identity, dict):
            return False
        path = Path(str(identity.get("path", "")))
        if file_identity(path) != identity:
            return False
    return True


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path(__file__).resolve().parents[3],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return ""


def parse_bool(value: object) -> Optional[bool]:
    if pd.isna(value):
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return None


def parse_float(value: object) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def evaluate_expression(
    row: pd.Series,
    min_tpm_tumor: float,
    max_tpm_ctrl: float,
    min_log2fc: float,
) -> tuple[bool, list[str], dict[str, object]]:
    tumor = parse_float(row.get("TPM_tumor", ""))
    ctrl = parse_float(row.get("TPM_ctrl", ""))
    log2fc = parse_float(row.get("log2FC", ""))
    if tumor is None or ctrl is None or log2fc is None:
        return (
            False,
            ["expression_value_missing_or_invalid"],
            {
                "TPM_tumor": "" if tumor is None else tumor,
                "TPM_ctrl": "" if ctrl is None else ctrl,
                "log2FC": "" if log2fc is None else log2fc,
                "expression_qc_status": "fail",
                "expression_qc_reasons": "expression_value_missing_or_invalid",
                "upstream_is_aberrant_consistent": "",
            },
        )
    reasons: list[str] = []
    if tumor < min_tpm_tumor:
        reasons.append("tumor_tpm_below_threshold")
    if ctrl > max_tpm_ctrl:
        reasons.append("control_tpm_above_threshold")
    if log2fc < min_log2fc:
        reasons.append("log2fc_below_threshold")
    passed = not reasons
    upstream = parse_bool(row.get("is_aberrant", ""))
    upstream_consistent = "" if upstream is None else bool(upstream == passed)
    if upstream is not None and upstream != passed:
        reasons.append("upstream_is_aberrant_mismatch")
    return (
        passed,
        [reason for reason in reasons if reason != "upstream_is_aberrant_mismatch"],
        {
            "TPM_tumor": tumor,
            "TPM_ctrl": ctrl,
            "log2FC": log2fc,
            "expression_qc_status": "pass" if passed else "fail",
            "expression_qc_reasons": ";".join(
                reason for reason in reasons if reason != "upstream_is_aberrant_mismatch"
            ),
            "upstream_is_aberrant": "" if upstream is None else upstream,
            "upstream_is_aberrant_consistent": upstream_consistent,
        },
    )


def stable_peptide_id(sample: str, mhc_class: str, peptide: str, index: int) -> str:
    digest = sha256_text(f"{sample}|{mhc_class}|{peptide}")[:16]
    class_tag = "MHC-I" if mhc_class == "MHC-I" else "MHC-II"
    return f"cryptic_core_{class_tag}_{index:07d}_{digest}"


def build_core(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    script_sha256 = sha256_file(Path(__file__))
    mhc_i_lengths = parse_int_list(args.mhc_i_lengths)
    mhc_ii_lengths = parse_int_list(args.mhc_ii_lengths)
    if set(mhc_i_lengths).intersection(mhc_ii_lengths):
        raise ValueError("MHC-I and MHC-II peptide lengths must not overlap for peptide-core binding")

    input_paths = {
        "aeseps_fasta": file_identity(Path(args.ae_seps_fasta)),
        "aeseps_annotation": file_identity(Path(args.aeseps_annotation)),
        "orf_filtered_fasta": file_identity(Path(args.orf_filtered_fasta)),
        "orf_final": file_identity(Path(args.orf_final)),
    }
    human_proteome = Path(args.human_proteome_fasta) if args.human_proteome_fasta else None
    if human_proteome:
        input_paths["human_reference_proteome_fasta"] = file_identity(human_proteome)
    elif not bool(args.allow_missing_human_reference):
        raise FileNotFoundError(
            "human reference proteome FASTA is required for formal Cryptic Core QC. "
            "Use --allow-missing-human-reference only for exploratory runs."
        )
    for key, value in [
        ("reference_genome_fasta", args.reference_genome_fasta),
        ("reference_gtf", args.reference_gtf),
        ("reference_lnc_gtf", args.reference_lnc_gtf),
    ]:
        if value:
            input_paths[key] = file_identity(Path(value))

    signature = {
        "policy_version": POLICY_VERSION,
        "script_sha256": script_sha256,
        "sample": args.sample,
        "matched_control_sample": args.matched_control_sample,
        "accepted_source_classes": sorted(ACCEPTED_SOURCE_CLASSES),
        "min_tpm_tumor": float(args.min_tpm_tumor),
        "max_tpm_ctrl": float(args.max_tpm_ctrl),
        "min_log2fc": float(args.min_log2fc),
        "mhc_i_lengths": mhc_i_lengths,
        "mhc_ii_lengths": mhc_ii_lengths,
        "allow_missing_human_reference": bool(args.allow_missing_human_reference),
        "strandedness": args.strandedness,
        "input_paths": input_paths,
    }

    manifest_path = outdir / "run_manifest.json"
    required_outputs = [
        outdir / "cryptic_parent_records.tsv",
        outdir / "cryptic_parent_core.tsv",
        outdir / "cryptic_parent_excluded.tsv",
        outdir / "cryptic_parent_core.fasta",
        outdir / "cryptic_peptide_core_hla_i.fasta",
        outdir / "cryptic_peptide_core_hla_ii.fasta",
        outdir / "cryptic_peptide_core.fasta",
        outdir / "cryptic_peptide_core_hla_i.tsv",
        outdir / "cryptic_peptide_core_hla_ii.tsv",
        outdir / "cryptic_peptide_core.tsv",
        outdir / "cryptic_peptide_parent_map.tsv",
        outdir / "cryptic_normal_resource_evidence.tsv",
        outdir / "stagewise_qc.tsv",
    ]
    if manifest_path.exists():
        old = read_json(manifest_path)
        if old.get("input_signature") != signature:
            raise ValueError(
                f"Existing Cryptic Core manifest does not match current inputs/config: {manifest_path}"
            )
        old_output_signature = old.get("output_signature")
        if isinstance(old_output_signature, dict) and outputs_match_signature(old_output_signature):
            return {**old, "reused": True}

    annot = pd.read_csv(args.aeseps_annotation, low_memory=False)
    if "Name" not in annot.columns or "nc_class" not in annot.columns or "is_aberrant" not in annot.columns:
        raise ValueError("aeSEP annotation table must contain Name, nc_class, and is_aberrant columns")
    annot["Name"] = annot["Name"].astype(str)
    annot_by_name = annot.drop_duplicates("Name", keep="first").set_index("Name")

    orf_final = pd.read_csv(args.orf_final)
    if "TranscriptID" not in orf_final.columns:
        raise ValueError("orf_final table must contain TranscriptID")
    allowed_orf_ids = set(orf_final["TranscriptID"].astype(str))

    parent_rows: list[dict[str, object]] = []
    parent_core_rows: list[dict[str, object]] = []
    parent_excluded_rows: list[dict[str, object]] = []
    candidate_peptide_rows: list[dict[str, object]] = []
    parent_fasta_rows: list[tuple[str, str]] = []

    seen_records: set[str] = set()
    for header, raw_seq in read_fasta(Path(args.orf_filtered_fasta)):
        record_id = fasta_record_id(header)
        if record_id in seen_records:
            raise ValueError(f"Duplicate ORF record id in ORF-filtered aeSEP FASTA: {record_id}")
        seen_records.add(record_id)
        src_id = source_parent_id(record_id)
        sequence, seq_reasons = normalize_parent_sequence(raw_seq)
        reasons = list(seq_reasons)

        if record_id not in allowed_orf_ids:
            reasons.append("not_in_orf_final")
        if src_id not in annot_by_name.index:
            reasons.append("source_join_failure")
            nc_class = ""
            expression_info = {
                "TPM_tumor": "",
                "TPM_ctrl": "",
                "log2FC": "",
                "expression_qc_status": "fail",
                "expression_qc_reasons": "source_join_failure",
                "upstream_is_aberrant": "",
                "upstream_is_aberrant_consistent": "",
            }
            tpm_tumor = ""
            tpm_ctrl = ""
            log2fc = ""
        else:
            row = annot_by_name.loc[src_id]
            nc_class = str(row.get("nc_class", ""))
            expression_pass, expression_reasons, expression_info = evaluate_expression(
                row,
                min_tpm_tumor=float(args.min_tpm_tumor),
                max_tpm_ctrl=float(args.max_tpm_ctrl),
                min_log2fc=float(args.min_log2fc),
            )
            tpm_tumor = expression_info["TPM_tumor"]
            tpm_ctrl = expression_info["TPM_ctrl"]
            log2fc = expression_info["log2FC"]
            if nc_class not in ACCEPTED_SOURCE_CLASSES:
                reasons.append(f"source_not_core:{nc_class or 'missing'}")
            reasons.extend(expression_reasons)

        parent_status = "core" if not reasons else "excluded"
        parent_row = {
            "sample": args.sample,
            "parent_record_id": record_id,
            "source_parent_id": src_id,
            "nc_class": nc_class,
            "upstream_is_aberrant": expression_info["upstream_is_aberrant"],
            "upstream_is_aberrant_consistent": expression_info["upstream_is_aberrant_consistent"],
            "TPM_tumor": tpm_tumor,
            "TPM_ctrl": tpm_ctrl,
            "log2FC": log2fc,
            "expression_qc_status": expression_info["expression_qc_status"],
            "expression_qc_reasons": expression_info["expression_qc_reasons"],
            "parent_sequence": sequence,
            "parent_sequence_sha256": sha256_text(sequence) if sequence else "",
            "parent_qc_status": parent_status,
            "parent_qc_reasons": ";".join(reasons),
            "rna_variant_status": "not_evaluated",
            "rna_editing_qc_status": "editing_resource_not_evaluated",
            "junction_support_status": "not_evaluated",
            "cryptic_discovery_fdr_status": "not_estimable_single_pair_rule_based_qc",
        }
        parent_rows.append(parent_row)
        if parent_status == "core":
            parent_core_rows.append(parent_row)
            parent_fasta_rows.append((record_id, sequence))
            for mhc_class, lengths in [("MHC-I", mhc_i_lengths), ("MHC-II", mhc_ii_lengths)]:
                for start, length, peptide in iter_windows(sequence, lengths):
                    candidate_peptide_rows.append({
                        "sample": args.sample,
                        "parent_record_id": record_id,
                        "source_parent_id": src_id,
                        "mhc_class": mhc_class,
                        "peptide_start": start,
                        "peptide_length": length,
                        "peptide": peptide,
                        "peptide_sequence_sha256": sha256_text(peptide),
                        "peptide_core_status": "core",
                        "peptide_qc_reasons": "",
                        "human_reference_peptide_status": "not_evaluated",
                        "normal_hla_presentation_status": "normal_hla_atlas_not_evaluated",
                        "normal_translation_status": "normal_translation_not_evaluable",
                        "external_normal_status": "external_normal_not_evaluated",
                    })
        else:
            parent_excluded_rows.append(parent_row)

    candidate_peptides = {str(row["peptide"]) for row in candidate_peptide_rows}
    human_matches, human_reference_summary = load_human_reference_matches(human_proteome, candidate_peptides)

    core_peptide_rows: list[dict[str, object]] = []
    excluded_peptide_rows: list[dict[str, object]] = []
    for row in candidate_peptide_rows:
        peptide = str(row["peptide"])
        if human_reference_summary["status"] == "evaluated":
            if peptide in human_matches:
                row["peptide_core_status"] = "excluded"
                row["peptide_qc_reasons"] = "human_reference_peptide_match"
                row["human_reference_peptide_status"] = "human_reference_peptide_match"
                excluded_peptide_rows.append(row)
                continue
            row["human_reference_peptide_status"] = "human_reference_peptide_not_detected"
        core_peptide_rows.append(row)

    peptide_id_by_key: dict[tuple[str, str], str] = {}
    unique_core_peptide_rows: list[dict[str, object]] = []
    support_counts: defaultdict[tuple[str, str], int] = defaultdict(int)
    for row in core_peptide_rows:
        key = (str(row["mhc_class"]), str(row["peptide"]))
        support_counts[key] += 1
        if key not in peptide_id_by_key:
            peptide_id = stable_peptide_id(args.sample, key[0], key[1], len(unique_core_peptide_rows) + 1)
            peptide_id_by_key[key] = peptide_id
            unique_row = dict(row)
            unique_row["peptide_record_id"] = peptide_id
            unique_row["supporting_window_count"] = 0
            unique_core_peptide_rows.append(unique_row)
        row["peptide_record_id"] = peptide_id_by_key[key]
    for row in unique_core_peptide_rows:
        key = (str(row["mhc_class"]), str(row["peptide"]))
        row["supporting_window_count"] = support_counts[key]
    for row in excluded_peptide_rows:
        row["peptide_record_id"] = ""
        row["supporting_window_count"] = ""

    hla_i_rows = [row for row in unique_core_peptide_rows if row["mhc_class"] == "MHC-I"]
    hla_ii_rows = [row for row in unique_core_peptide_rows if row["mhc_class"] == "MHC-II"]

    with (outdir / "cryptic_parent_core.fasta").open("w") as handle:
        for record_id, sequence in parent_fasta_rows:
            handle.write(f">{record_id}\n{sequence}\n")
    write_peptide_core_fasta(hla_i_rows, outdir / "cryptic_peptide_core_hla_i.fasta")
    write_peptide_core_fasta(hla_ii_rows, outdir / "cryptic_peptide_core_hla_ii.fasta")
    write_peptide_core_fasta(unique_core_peptide_rows, outdir / "cryptic_peptide_core.fasta")

    parent_fields = [
        "sample", "parent_record_id", "source_parent_id", "nc_class",
        "upstream_is_aberrant", "upstream_is_aberrant_consistent",
        "TPM_tumor", "TPM_ctrl", "log2FC", "expression_qc_status", "expression_qc_reasons",
        "parent_sequence", "parent_sequence_sha256",
        "parent_qc_status", "parent_qc_reasons", "rna_variant_status",
        "rna_editing_qc_status", "junction_support_status", "cryptic_discovery_fdr_status",
    ]
    peptide_fields = [
        "sample", "peptide_record_id", "parent_record_id", "source_parent_id", "mhc_class", "peptide_start",
        "peptide_length", "peptide", "peptide_sequence_sha256", "peptide_core_status",
        "peptide_qc_reasons", "human_reference_peptide_status",
        "normal_hla_presentation_status", "normal_translation_status", "external_normal_status",
        "supporting_window_count",
    ]
    write_tsv(parent_rows, outdir / "cryptic_parent_records.tsv", parent_fields)
    write_tsv(parent_core_rows, outdir / "cryptic_parent_core.tsv", parent_fields)
    write_tsv(parent_excluded_rows, outdir / "cryptic_parent_excluded.tsv", parent_fields)
    write_tsv(hla_i_rows, outdir / "cryptic_peptide_core_hla_i.tsv", peptide_fields)
    write_tsv(hla_ii_rows, outdir / "cryptic_peptide_core_hla_ii.tsv", peptide_fields)
    write_tsv(unique_core_peptide_rows, outdir / "cryptic_peptide_core.tsv", peptide_fields)
    write_tsv(core_peptide_rows + excluded_peptide_rows, outdir / "cryptic_peptide_parent_map.tsv", peptide_fields)
    write_tsv(core_peptide_rows + excluded_peptide_rows, outdir / "cryptic_normal_resource_evidence.tsv", peptide_fields)

    stage_counts = [
        {"stage": "orf_filtered_parent_records", "count": len(parent_rows)},
        {"stage": "parent_core_records", "count": len(parent_core_rows)},
        {"stage": "parent_excluded_records", "count": len(parent_excluded_rows)},
        {"stage": "candidate_peptide_window_rows", "count": len(candidate_peptide_rows)},
        {"stage": "peptide_core_window_rows", "count": len(core_peptide_rows)},
        {"stage": "peptide_excluded_window_rows", "count": len(excluded_peptide_rows)},
        {"stage": "unique_peptide_core_records", "count": len(unique_core_peptide_rows)},
        {"stage": "hla_i_unique_peptide_core_records", "count": len(hla_i_rows)},
        {"stage": "hla_ii_unique_peptide_core_records", "count": len(hla_ii_rows)},
    ]
    write_tsv(stage_counts, outdir / "stagewise_qc.tsv", ["stage", "count"])

    outputs_hash = output_signature(required_outputs)
    manifest = {
        "sample": args.sample,
        "matched_control_sample": args.matched_control_sample,
        "policy_version": POLICY_VERSION,
        "run_status": "complete",
        "started_at": ts(),
        "finished_at": ts(),
        "code_commit": git_commit(),
        "script_sha256": script_sha256,
        "input_signature": signature,
        "human_reference_evaluation": human_reference_summary["status"],
        "human_reference_summary": human_reference_summary,
        "outputs": {path.name: str(path) for path in required_outputs},
        "output_signature": outputs_hash,
        "stage_counts": {row["stage"]: row["count"] for row in stage_counts},
        "binding_input_fasta": str(outdir / "cryptic_peptide_core.fasta"),
        "binding_input_mode": "peptide-core",
        "reused": False,
    }
    write_json(manifest_path, manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", required=True)
    parser.add_argument("--matched-control-sample", default="")
    parser.add_argument("--ae-seps-fasta", required=True)
    parser.add_argument("--aeseps-annotation", required=True)
    parser.add_argument("--orf-filtered-fasta", required=True)
    parser.add_argument("--orf-final", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--human-proteome-fasta", default="")
    parser.add_argument(
        "--allow-missing-human-reference",
        action="store_true",
        help="Allow exploratory Core generation without human reference exact peptide matching.",
    )
    parser.add_argument("--reference-genome-fasta", default="")
    parser.add_argument("--reference-gtf", default="")
    parser.add_argument("--reference-lnc-gtf", default="")
    parser.add_argument("--strandedness", default="")
    parser.add_argument("--min-tpm-tumor", type=float, default=5.0)
    parser.add_argument("--max-tpm-ctrl", type=float, default=0.5)
    parser.add_argument("--min-log2fc", type=float, default=4.0)
    parser.add_argument("--mhc-i-lengths", default="8,9,10,11")
    parser.add_argument("--mhc-ii-lengths", default="13,14,15,16,17")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    manifest = build_core(build_parser().parse_args(argv))
    print(json.dumps(manifest, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
