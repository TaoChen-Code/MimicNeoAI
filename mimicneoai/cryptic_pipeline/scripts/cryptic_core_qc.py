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

try:
    from mimicneoai.cryptic_pipeline.scripts.cryptic_coordinate_utils import (
        GenomicBlock,
        blocks_to_json,
        blocks_to_text,
        canonical_chromosome,
        extract_transcript_cds_from_reference,
        genomic_order,
        load_bam_alignment_stats,
        map_transcript_interval_to_genome,
        open_reference_fasta,
        parse_bed12,
        reference_contig_lookup,
        read_fasta_records,
        transcript_order,
        translate_cds,
    )
except ImportError:  # pragma: no cover - direct script execution fallback
    from cryptic_coordinate_utils import (  # type: ignore
        GenomicBlock,
        blocks_to_json,
        blocks_to_text,
        canonical_chromosome,
        extract_transcript_cds_from_reference,
        genomic_order,
        load_bam_alignment_stats,
        map_transcript_interval_to_genome,
        open_reference_fasta,
        parse_bed12,
        reference_contig_lookup,
        read_fasta_records,
        transcript_order,
        translate_cds,
    )

POLICY_VERSION_V1 = "cryptic_core_qc_v1.0"
POLICY_VERSION_V11 = "cryptic_core_qc_v1.1"
POLICY_VERSION = POLICY_VERSION_V1
COORDINATE_UTILS_PATH = Path(__file__).with_name("cryptic_coordinate_utils.py")
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


def iter_canonical_segments(sequence: str) -> Iterable[str]:
    start: Optional[int] = None
    for index, aa in enumerate(sequence):
        if aa in CANONICAL_AA:
            if start is None:
                start = index
        elif start is not None:
            yield sequence[start:index]
            start = None
    if start is not None:
        yield sequence[start:]


def normalize_human_reference_sequence(sequence: str) -> str:
    seq = str(sequence).replace(" ", "").replace("\n", "").strip().upper()
    if seq.endswith("*"):
        seq = seq[:-1]
    return seq


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
    records_encountered = 0
    records_with_noncanonical_residues = 0
    standard_windows_evaluated = 0
    for _header, raw_seq in read_fasta(human_proteome_fasta):
        records_encountered += 1
        seq = normalize_human_reference_sequence(raw_seq)
        if not seq:
            continue
        if set(seq).difference(CANONICAL_AA):
            records_with_noncanonical_residues += 1
        for segment in iter_canonical_segments(seq):
            for length in lengths:
                if len(segment) < length:
                    continue
                standard_windows_evaluated += len(segment) - length + 1
                for start in range(0, len(segment) - length + 1):
                    pep = segment[start : start + length]
                    if pep in candidate_peptides:
                        matches.add(pep)
    return matches, {
        "status": "evaluated",
        "path": str(human_proteome_fasta),
        "sha256": sha256_file(human_proteome_fasta),
        "records_scanned": records_encountered,
        "records_encountered": records_encountered,
        "records_with_noncanonical_residues": records_with_noncanonical_residues,
        "standard_windows_evaluated": standard_windows_evaluated,
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


def is_v11(policy_version: str) -> bool:
    return policy_version == POLICY_VERSION_V11


def require_coordinate_inputs(args: argparse.Namespace) -> None:
    required = [
        ("--orf-bed12", args.orf_bed12),
        ("--orf-bam", args.orf_bam),
        ("--orf-cds-fasta", args.orf_cds_fasta),
        ("--reference-genome-fasta", args.reference_genome_fasta),
    ]
    missing = [label for label, value in required if not str(value or "").strip()]
    if missing:
        raise ValueError(
            "cryptic_core_qc_v1.1 requires coordinate inputs: " + ", ".join(missing)
        )
    for label, value in required:
        path = Path(value)
        if not path.exists():
            raise FileNotFoundError(f"{label} file not found: {path}")


def summarize_alignment_loci(alignments: list[dict[str, object]]) -> str:
    values = []
    for aln in alignments:
        chrom = str(aln.get("chromosome", ""))
        strand = str(aln.get("strand", ""))
        start = str(aln.get("start0", ""))
        end = str(aln.get("end0", ""))
        mapq = str(aln.get("mapq", ""))
        cigar = str(aln.get("cigar", ""))
        values.append(f"{chrom}:{start}-{end}:{strand}:MAPQ{mapq}:{cigar}")
    return ";".join(values)


def build_coordinate_sidecars(
    args: argparse.Namespace,
    parent_core_rows: list[dict[str, object]],
    peptide_parent_rows: list[dict[str, object]],
    policy_version: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    if not is_v11(policy_version):
        return [], [], [], {"status": "not_evaluated_policy_v1.0"}

    require_coordinate_inputs(args)
    bed12_entries = parse_bed12(Path(args.orf_bed12))
    bam_stats = load_bam_alignment_stats(Path(args.orf_bam))
    cds_records = read_fasta_records(Path(args.orf_cds_fasta))
    reference_fasta = open_reference_fasta(Path(args.reference_genome_fasta))
    reference_lookup = reference_contig_lookup(reference_fasta)
    min_mapq = int(args.coordinate_min_mapq)
    reference_build = str(args.reference_build or "GRCh38")

    parent_by_id = {str(row["parent_record_id"]): row for row in parent_core_rows}
    parent_coordinates: list[dict[str, object]] = []
    parent_orfcds: list[dict[str, object]] = []
    parent_status: dict[str, dict[str, object]] = {}

    try:
        for parent_id, parent_row in parent_by_id.items():
            stats = bam_stats.get(parent_id, {"primary": [], "secondary": [], "supplementary": []})
            primary = list(stats.get("primary", []))
            secondary = list(stats.get("secondary", []))
            supplementary = list(stats.get("supplementary", []))
            bed_entries = bed12_entries.get(parent_id, [])
            reasons: list[str] = []
            status = "coordinate_evaluable"
            chosen: dict[str, object] | None = primary[0] if len(primary) == 1 else None

            if len(primary) != 1:
                status = "ambiguous_candidate_mapping"
                reasons.append("primary_alignment_count_not_one")
            if secondary:
                status = "ambiguous_candidate_mapping"
                reasons.append("secondary_alignment_present")
            if supplementary:
                status = "ambiguous_candidate_mapping"
                reasons.append("supplementary_alignment_present")
            if chosen is None:
                chosen = primary[0] if primary else None

            blocks = [block for block in (chosen or {}).get("blocks", []) if isinstance(block, GenomicBlock)]
            chrom = str((chosen or {}).get("chromosome", ""))
            strand = str((chosen or {}).get("strand", ""))
            mapq_values = [int(aln.get("mapq", 0)) for aln in primary]
            mapq = mapq_values[0] if len(mapq_values) == 1 else ""
            if chosen is not None and mapq != "" and int(mapq) < min_mapq:
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.append(f"mapq_below_{min_mapq}")
            if chosen is not None and not canonical_chromosome(chrom):
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.append("noncanonical_contig")

            tx_blocks = transcript_order(blocks, strand) if blocks and strand in {"+", "-"} else []
            g_blocks = genomic_order(blocks)
            input_cds_seq = cds_records.get(parent_id, "")
            input_cds_length = len(input_cds_seq)
            block_total = sum(block.length for block in blocks)
            input_translated = translate_cds(input_cds_seq) if input_cds_seq else ""
            parent_sequence = str(parent_row.get("parent_sequence", ""))
            input_translation_status = "pass"
            if not input_cds_seq:
                input_translation_status = "missing_cds_fasta_record"
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.append("missing_cds_fasta_record")
            elif input_translated != parent_sequence:
                input_translation_status = "input_cds_translation_mismatch"
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.append("input_cds_translation_mismatch")
            if chosen is not None and input_cds_seq and block_total != input_cds_length:
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.append("input_cds_length_block_total_mismatch")

            reference_cds_seq = ""
            reference_reasons: list[str] = []
            reference_translation_status = "not_evaluable_candidate_coordinate_qc"
            rna_variant_coordinate_status = "not_evaluated"
            if chosen is None or not tx_blocks:
                reference_reasons.append("missing_trusted_primary_blocks")
            else:
                reference_cds_seq, reference_reasons = extract_transcript_cds_from_reference(
                    reference_fasta,
                    reference_lookup,
                    tx_blocks,
                )
            if reference_reasons:
                status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                reasons.extend(reference_reasons)
                reference_translation_status = "reference_translation_not_evaluable"
            else:
                reference_translated = translate_cds(reference_cds_seq)
                if reference_translated == parent_sequence:
                    reference_translation_status = "pass"
                else:
                    status = "not_evaluable_candidate_coordinate_qc" if status == "coordinate_evaluable" else status
                    reasons.append("reference_translation_mismatch")
                    reference_translation_status = "reference_translation_mismatch"
                    rna_variant_coordinate_status = "RNA_variant_aware_not_evaluated"
            terminal_stop = (
                len(input_cds_seq) >= 3
                and len(input_cds_seq) % 3 == 0
                and input_cds_seq[-3:].upper() in {"TAA", "TAG", "TGA"}
            )
            parent_info = {
                "sample": args.sample,
                "parent_record_id": parent_id,
                "chromosome": chrom,
                "strand": strand,
                "genomic_start0": min((block.start0 for block in blocks), default=""),
                "genomic_end0": max((block.end0 for block in blocks), default=""),
                "block_count": len(blocks),
                "source_blocks_genomic_order": blocks_to_text(g_blocks),
                "blocks_transcript_order": blocks_to_text(tx_blocks),
                "primary_alignment_count": len(primary),
                "secondary_alignment_count": len(secondary),
                "supplementary_alignment_count": len(supplementary),
                "MAPQ": mapq,
                "cds_nucleotide_length": input_cds_length,
                "block_total_length": block_total,
                "terminal_stop_in_cds": terminal_stop,
                "reference_genomic_cds_length": len(reference_cds_seq),
                "input_cds_translation_status": input_translation_status,
                "reference_genomic_translation_status": reference_translation_status,
                "rna_variant_coordinate_status": rna_variant_coordinate_status,
                "reference_build": reference_build,
                "bed12_alignment_count": len(bed_entries),
                "all_alignment_loci": summarize_alignment_loci(primary + secondary + supplementary),
                "coordinate_mapping_status": status,
                "coordinate_mapping_reasons": ";".join(sorted(set(reasons))),
            }
            parent_coordinates.append(parent_info)
            parent_status[parent_id] = {**parent_info, "transcript_blocks": tx_blocks}

            for idx, block in enumerate(tx_blocks, start=1):
                cumulative_before = sum(prev.length for prev in tx_blocks[: idx - 1])
                parent_orfcds.append({
                    "sample": args.sample,
                    "parent_record_id": parent_id,
                    "transcript_block_order": idx,
                    "chromosome": block.chrom,
                    "strand": block.strand,
                    "start0": block.start0,
                    "end0": block.end0,
                    "block_length": block.length,
                    "derived_phase": cumulative_before % 3,
                    "phase_provenance": "derived_from_zero_phase_candidate_orf_start",
                    "coordinate_mapping_status": status,
                })
    finally:
        reference_fasta.close()

    peptide_footprints: list[dict[str, object]] = []
    for row in peptide_parent_rows:
        parent_id = str(row.get("parent_record_id", ""))
        parent_info = parent_status.get(parent_id, {})
        status = str(parent_info.get("coordinate_mapping_status", "not_evaluable_candidate_coordinate_qc"))
        reasons = str(parent_info.get("coordinate_mapping_reasons", ""))
        tx_blocks = parent_info.get("transcript_blocks", [])
        peptide_start = int(row.get("peptide_start", 0))
        peptide_length = int(row.get("peptide_length", 0))
        cds_start0 = (peptide_start - 1) * 3
        cds_end0 = cds_start0 + peptide_length * 3
        footprint_blocks: list[GenomicBlock] = []
        if status == "coordinate_evaluable":
            try:
                footprint_blocks = [
                    block
                    for _segment_start, _segment_end, block in map_transcript_interval_to_genome(
                        tx_blocks, str(parent_info.get("strand", "")), cds_start0, cds_end0
                    )
                ]
            except Exception as exc:
                status = "not_evaluable_candidate_coordinate_qc"
                reasons = ";".join(filter(None, [reasons, f"peptide_footprint_mapping_failed:{exc}"]))
        peptide_footprints.append({
            "sample": args.sample,
            "peptide_record_id": row.get("peptide_record_id", ""),
            "parent_record_id": parent_id,
            "source_parent_id": row.get("source_parent_id", ""),
            "mhc_class": row.get("mhc_class", ""),
            "peptide": row.get("peptide", ""),
            "peptide_start": peptide_start,
            "peptide_length": peptide_length,
            "chromosome": parent_info.get("chromosome", ""),
            "strand": parent_info.get("strand", ""),
            "peptide_cds_start0": cds_start0,
            "peptide_cds_end0": cds_end0,
            "codon_blocks_transcript_order": blocks_to_json(footprint_blocks),
            "codon_blocks_text": blocks_to_text(footprint_blocks),
            "peptide_footprint_phase": cds_start0 % 3,
            "candidate_coordinate_status": status,
            "candidate_coordinate_reasons": reasons,
            "reference_build": reference_build,
        })

    summary = {
        "status": "evaluated",
        "parent_rows": len(parent_coordinates),
        "parent_coordinate_evaluable": sum(
            1 for row in parent_coordinates if row["coordinate_mapping_status"] == "coordinate_evaluable"
        ),
        "parent_coordinate_not_evaluable": sum(
            1 for row in parent_coordinates if row["coordinate_mapping_status"] != "coordinate_evaluable"
        ),
        "parent_reference_translation_pass": sum(
            1 for row in parent_coordinates if row["reference_genomic_translation_status"] == "pass"
        ),
        "parent_reference_translation_mismatch": sum(
            1 for row in parent_coordinates
            if row["reference_genomic_translation_status"] == "reference_translation_mismatch"
        ),
        "parent_reference_translation_not_evaluable": sum(
            1 for row in parent_coordinates
            if row["reference_genomic_translation_status"] == "reference_translation_not_evaluable"
        ),
        "peptide_footprint_rows": len(peptide_footprints),
        "peptide_coordinate_evaluable": sum(
            1 for row in peptide_footprints if row["candidate_coordinate_status"] == "coordinate_evaluable"
        ),
    }
    return parent_coordinates, parent_orfcds, peptide_footprints, summary


def build_core(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    policy_version = str(args.policy_version)
    if policy_version not in {POLICY_VERSION_V1, POLICY_VERSION_V11}:
        raise ValueError(f"Unsupported Cryptic Core QC policy_version: {policy_version}")
    if is_v11(policy_version):
        require_coordinate_inputs(args)

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
    if is_v11(policy_version):
        input_paths["orf_bed12"] = file_identity(Path(args.orf_bed12))
        input_paths["orf_bam"] = file_identity(Path(args.orf_bam))
        input_paths["orf_cds_fasta"] = file_identity(Path(args.orf_cds_fasta))
        input_paths["coordinate_utils"] = file_identity(COORDINATE_UTILS_PATH)
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
        "policy_version": policy_version,
        "script_sha256": script_sha256,
        "coordinate_utils_sha256": sha256_file(COORDINATE_UTILS_PATH) if is_v11(policy_version) else "",
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
        "coordinate_min_mapq": int(args.coordinate_min_mapq),
        "reference_build": str(args.reference_build or "GRCh38"),
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
    if is_v11(policy_version):
        required_outputs.extend([
            outdir / "cryptic_parent_coordinates.tsv",
            outdir / "cryptic_parent_orfcds.tsv",
            outdir / "cryptic_peptide_genomic_footprint.tsv",
        ])
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
    parent_coordinate_rows, parent_orfcds_rows, peptide_footprint_rows, coordinate_summary = build_coordinate_sidecars(
        args,
        parent_core_rows,
        core_peptide_rows,
        policy_version,
    )

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
    if is_v11(policy_version):
        parent_coordinate_fields = [
            "sample", "parent_record_id", "chromosome", "strand", "genomic_start0", "genomic_end0",
            "block_count", "source_blocks_genomic_order", "blocks_transcript_order",
            "primary_alignment_count", "secondary_alignment_count", "supplementary_alignment_count",
            "MAPQ", "cds_nucleotide_length", "block_total_length", "terminal_stop_in_cds",
            "reference_genomic_cds_length", "input_cds_translation_status",
            "reference_genomic_translation_status", "rna_variant_coordinate_status",
            "reference_build", "bed12_alignment_count", "all_alignment_loci",
            "coordinate_mapping_status", "coordinate_mapping_reasons",
        ]
        parent_orfcds_fields = [
            "sample", "parent_record_id", "transcript_block_order", "chromosome", "strand",
            "start0", "end0", "block_length", "derived_phase", "phase_provenance",
            "coordinate_mapping_status",
        ]
        peptide_footprint_fields = [
            "sample", "peptide_record_id", "parent_record_id", "source_parent_id", "mhc_class",
            "peptide", "peptide_start", "peptide_length", "chromosome", "strand",
            "peptide_cds_start0", "peptide_cds_end0", "codon_blocks_transcript_order",
            "codon_blocks_text", "peptide_footprint_phase", "candidate_coordinate_status",
            "candidate_coordinate_reasons", "reference_build",
        ]
        write_tsv(parent_coordinate_rows, outdir / "cryptic_parent_coordinates.tsv", parent_coordinate_fields)
        write_tsv(parent_orfcds_rows, outdir / "cryptic_parent_orfcds.tsv", parent_orfcds_fields)
        write_tsv(peptide_footprint_rows, outdir / "cryptic_peptide_genomic_footprint.tsv", peptide_footprint_fields)

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
    if is_v11(policy_version):
        stage_counts.extend([
            {
                "stage": "candidate_parents_coordinate_evaluable",
                "count": int(coordinate_summary.get("parent_coordinate_evaluable", 0)),
            },
            {
                "stage": "candidate_parents_coordinate_not_evaluable",
                "count": int(coordinate_summary.get("parent_coordinate_not_evaluable", 0)),
            },
            {
                "stage": "candidate_parents_reference_translation_pass",
                "count": int(coordinate_summary.get("parent_reference_translation_pass", 0)),
            },
            {
                "stage": "candidate_parents_reference_translation_mismatch",
                "count": int(coordinate_summary.get("parent_reference_translation_mismatch", 0)),
            },
            {
                "stage": "candidate_parents_reference_translation_not_evaluable",
                "count": int(coordinate_summary.get("parent_reference_translation_not_evaluable", 0)),
            },
            {
                "stage": "peptides_coordinate_evaluable",
                "count": int(coordinate_summary.get("peptide_coordinate_evaluable", 0)),
            },
        ])
    write_tsv(stage_counts, outdir / "stagewise_qc.tsv", ["stage", "count"])

    outputs_hash = output_signature(required_outputs)
    manifest = {
        "sample": args.sample,
        "matched_control_sample": args.matched_control_sample,
        "policy_version": policy_version,
        "run_status": "complete",
        "started_at": ts(),
        "finished_at": ts(),
        "code_commit": git_commit(),
        "script_sha256": script_sha256,
        "input_signature": signature,
        "human_reference_evaluation": human_reference_summary["status"],
        "human_reference_summary": human_reference_summary,
        "candidate_coordinate_summary": coordinate_summary,
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
    parser.add_argument(
        "--policy-version",
        default=POLICY_VERSION_V1,
        choices=[POLICY_VERSION_V1, POLICY_VERSION_V11],
        help="Cryptic Core QC policy version. v1.1 writes candidate coordinate sidecars.",
    )
    parser.add_argument("--matched-control-sample", default="")
    parser.add_argument("--ae-seps-fasta", required=True)
    parser.add_argument("--aeseps-annotation", required=True)
    parser.add_argument("--orf-filtered-fasta", required=True)
    parser.add_argument("--orf-final", required=True)
    parser.add_argument("--orf-bed12", default="")
    parser.add_argument("--orf-bam", default="")
    parser.add_argument("--orf-cds-fasta", default="")
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
    parser.add_argument("--reference-build", default="GRCh38")
    parser.add_argument("--strandedness", default="")
    parser.add_argument("--min-tpm-tumor", type=float, default=5.0)
    parser.add_argument("--max-tpm-ctrl", type=float, default=0.5)
    parser.add_argument("--min-log2fc", type=float, default=4.0)
    parser.add_argument("--mhc-i-lengths", default="8,9,10,11")
    parser.add_argument("--mhc-ii-lengths", default="13,14,15,16,17")
    parser.add_argument("--coordinate-min-mapq", type=int, default=20)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    manifest = build_core(build_parser().parse_args(argv))
    print(json.dumps(manifest, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
