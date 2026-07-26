#!/usr/bin/env python3
"""Build pair-level microbial peptide Core after matched-normal subtraction.

This step consumes per-sample protein-hit tables from
06.MicrobialPeptidesIdentification and writes a tumor-only peptide Core for
binding prediction. Matched-normal subtraction uses all technically valid normal
protein hits; contaminant blacklist filtering is applied afterward to tumor
residual peptides. It does not run HLA typing or binding prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from mimicneoai.microbial_pipeline.scripts.get_data_for_binding_pred import (
    extract_qseqid_tag,
    normalize_parent_sequence,
)


POLICY_VERSION = "microbial_pair_core_qc_v1.0"
STANDARD_COLUMNS = {
    "qseqid",
    "sseqid",
    "sseq",
    "pident",
    "evalue",
    "qcovs",
}
MHC_CLASS_I = "MHC-I"
MHC_CLASS_II = "MHC-II"


def parse_int_list(value: str) -> list[int]:
    lengths = []
    for token in str(value).replace(",", " ").split():
        token = token.strip()
        if token:
            lengths.append(int(token))
    if not lengths:
        raise ValueError("At least one peptide length is required")
    return sorted(set(lengths))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_sseqid(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    if text.startswith("ref|"):
        text = text[4:]
    return text.rstrip("|")


def split_taxids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [token.strip() for token in text.split(",") if token.strip()]


def validate_pair_id(pair_id: str, tumor_sample: str, normal_sample: str) -> None:
    if not tumor_sample or not normal_sample:
        raise ValueError("tumor_sample and normal_sample must be non-empty")
    if tumor_sample == normal_sample:
        raise ValueError("tumor_sample and normal_sample must be different")
    if pair_id:
        parts = [part.strip() for part in pair_id.split(",")]
        if len(parts) != 2 or any(not part for part in parts):
            raise ValueError("pair_id must be formatted as Tumor,Normal")
        if parts != [tumor_sample, normal_sample]:
            raise ValueError("pair_id does not match tumor_sample and normal_sample")


def load_contaminant_taxids(
    path: Optional[Path],
    *,
    allow_missing: bool = False,
    expected_sha256: str = "",
) -> tuple[set[str], dict[str, str], str]:
    if path is None:
        if allow_missing:
            return set(), {}, ""
        raise FileNotFoundError(
            "Contaminant blacklist is required for paired microbial Core QC. "
            "Use --allow-missing-blacklist only for exploratory runs."
        )
    if not path.exists():
        raise FileNotFoundError(f"Contaminant blacklist not found: {path}")
    actual_sha256 = sha256_file(path)
    if expected_sha256 and actual_sha256 != expected_sha256:
        raise ValueError(
            f"Contaminant blacklist SHA256 mismatch for {path}: "
            f"expected {expected_sha256}, observed {actual_sha256}"
        )
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"tax_id", "verdict"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Contaminant blacklist missing required columns: {missing}")
    verdict_by_taxid = {
        str(row.tax_id).strip(): str(row.verdict).strip()
        for row in df.itertuples(index=False)
        if str(row.tax_id).strip()
    }
    contaminants = {
        taxid
        for taxid, verdict in verdict_by_taxid.items()
        if verdict.lower() == "contaminant"
    }
    return contaminants, verdict_by_taxid, actual_sha256


def classify_blacklist_status(taxids: Iterable[str], contaminant_taxids: set[str]) -> str:
    taxid_set = {str(t).strip() for t in taxids if str(t).strip()}
    if not contaminant_taxids:
        return "not_evaluated"
    if not taxid_set:
        return "unresolved"
    contaminant_overlap = taxid_set.intersection(contaminant_taxids)
    if not contaminant_overlap:
        return "not_contaminant"
    if contaminant_overlap == taxid_set:
        return "exclusive_contaminant"
    return "mixed"


def read_protein_hits(
    path: Path,
    sample: str,
    role: str,
    min_pident: float,
    max_evalue: float,
    min_qcovs: float,
    contaminant_taxids: set[str],
    exclude_blacklist: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read new or legacy protein-hit table and split pass/excluded rows."""
    if not path.exists():
        raise FileNotFoundError(f"{role} protein-hit table not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        columns = [
            "parent_record_id",
            "sample",
            "role",
            "source_row_number",
            "qseqid",
            "protein_accession",
            "ctaxa_ids",
            "canonical_parent_sequence",
            "parent_sequence_sha256",
            "pident",
            "evalue",
            "qcovs",
            "blacklist_status",
            "parent_qc_status",
            "parent_qc_reasons",
            "parent_qc_flags",
        ]
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=columns)

    missing_standard = sorted(STANDARD_COLUMNS.difference(df.columns))
    if missing_standard and "canonical_parent_sequence" not in df.columns:
        raise ValueError(
            f"{role} protein-hit table is neither standard nor legacy-compatible; "
            f"missing columns: {missing_standard}"
        )
    if "qcovs" not in df.columns:
        raise ValueError(f"{role} protein-hit table missing required column: qcovs")

    source_numbers = (
        pd.to_numeric(df.get("source_row_number", pd.Series(range(1, len(df) + 1))), errors="coerce")
        .fillna(pd.Series(range(1, len(df) + 1)))
        .astype(int)
    )

    if "canonical_parent_sequence" in df.columns:
        raw_sequence = df["canonical_parent_sequence"]
    else:
        raw_sequence = df["sseq"]

    normalized = raw_sequence.map(normalize_parent_sequence)
    canonical_sequence = [item[0] for item in normalized]
    sequence_reasons = [item[1] for item in normalized]

    pident = pd.to_numeric(df.get("pident"), errors="coerce")
    evalue = pd.to_numeric(df.get("evalue"), errors="coerce")
    qcovs = pd.to_numeric(df.get("qcovs"), errors="coerce")

    protein_accession = df.get("protein_accession", pd.Series([""] * len(df))).fillna("").astype(str)
    if protein_accession.eq("").all() and "qseqid" in df.columns:
        protein_accession = df["qseqid"].map(lambda value: extract_qseqid_tag(value, "PROT"))
    if protein_accession.eq("").all() and "sseqid" in df.columns:
        protein_accession = df["sseqid"].map(normalize_sseqid)

    ctaxa_ids = df.get("ctaxa_ids", pd.Series([""] * len(df))).fillna("").astype(str)
    if ctaxa_ids.eq("").all() and "qseqid" in df.columns:
        ctaxa_ids = df["qseqid"].map(lambda value: extract_qseqid_tag(value, "CTAXA"))

    rows = []
    for i in range(len(df)):
        reasons = list(sequence_reasons[i])
        if pd.isna(pident.iloc[i]) or float(pident.iloc[i]) != float(min_pident):
            reasons.append(f"pident_not_equal_{min_pident:g}")
        if pd.isna(evalue.iloc[i]) or float(evalue.iloc[i]) > float(max_evalue):
            reasons.append(f"evalue_gt_{max_evalue:g}")
        if pd.isna(qcovs.iloc[i]) or float(qcovs.iloc[i]) < float(min_qcovs):
            reasons.append(f"qcovs_lt_{min_qcovs:g}")

        taxids = split_taxids(ctaxa_ids.iloc[i])
        blacklist_status = classify_blacklist_status(taxids, contaminant_taxids)
        flags = []
        if blacklist_status == "exclusive_contaminant":
            if exclude_blacklist:
                reasons.append("exclusive_contaminant")
            else:
                flags.append("blacklist_exclusive_contaminant_retained_for_subtraction")
        elif blacklist_status in {"mixed", "unresolved"}:
            flags.append(f"blacklist_{blacklist_status}_retained")

        source_row_number = int(source_numbers.iloc[i])
        sequence = canonical_sequence[i]
        parent_hash = sha256_text(sequence) if sequence else ""
        parent_record_id = f"{role}:{sample}:{source_row_number}"
        row = {
            "parent_record_id": parent_record_id,
            "sample": sample,
            "role": role,
            "source_row_number": source_row_number,
            "qseqid": df["qseqid"].iloc[i] if "qseqid" in df.columns else "",
            "protein_accession": protein_accession.iloc[i],
            "ctaxa_ids": ctaxa_ids.iloc[i],
            "canonical_parent_sequence": sequence,
            "parent_sequence_sha256": parent_hash,
            "pident": "" if pd.isna(pident.iloc[i]) else float(pident.iloc[i]),
            "evalue": "" if pd.isna(evalue.iloc[i]) else float(evalue.iloc[i]),
            "qcovs": "" if pd.isna(qcovs.iloc[i]) else float(qcovs.iloc[i]),
            "blacklist_status": blacklist_status,
            "parent_qc_status": "excluded" if reasons else "pass",
            "parent_qc_reasons": ";".join(reasons),
            "parent_qc_flags": ";".join(flags),
        }
        rows.append(row)

    all_rows = pd.DataFrame(rows)
    passing = all_rows[all_rows["parent_qc_status"] == "pass"].copy()
    excluded = all_rows[all_rows["parent_qc_status"] != "pass"].copy()
    return passing.reset_index(drop=True), excluded.reset_index(drop=True)


def iter_peptides(sequence: str, lengths: Iterable[int]) -> Iterable[tuple[str, int]]:
    for length in lengths:
        if len(sequence) < length:
            continue
        for start in range(0, len(sequence) - length + 1):
            yield sequence[start : start + length], start + 1


def build_peptide_maps(
    parents: pd.DataFrame,
    mhc_i_lengths: list[int],
    mhc_ii_lengths: list[int],
) -> tuple[dict[tuple[str, str], dict], dict[tuple[str, str, str], list[int]]]:
    """Build unique peptide map and peptide-parent position map."""
    peptide_map: dict[tuple[str, str], dict] = {}
    parent_positions: dict[tuple[str, str, str], list[int]] = defaultdict(list)

    for row in parents.itertuples(index=False):
        parent_id = str(row.parent_record_id)
        sequence = str(row.canonical_parent_sequence)
        for mhc_class, lengths in ((MHC_CLASS_I, mhc_i_lengths), (MHC_CLASS_II, mhc_ii_lengths)):
            for peptide, start in iter_peptides(sequence, lengths):
                key = (mhc_class, peptide)
                entry = peptide_map.setdefault(
                    key,
                    {
                        "mhc_class": mhc_class,
                        "peptide": peptide,
                        "peptide_length": len(peptide),
                        "parent_record_ids": set(),
                        "protein_accessions": set(),
                        "ctaxa_ids": set(),
                        "blacklist_statuses": set(),
                    },
                )
                entry["parent_record_ids"].add(parent_id)
                if str(row.protein_accession):
                    entry["protein_accessions"].add(str(row.protein_accession))
                for taxid in split_taxids(row.ctaxa_ids):
                    entry["ctaxa_ids"].add(taxid)
                entry["blacklist_statuses"].add(str(row.blacklist_status))
                parent_positions[(mhc_class, peptide, parent_id)].append(start)

    return peptide_map, parent_positions


def peptide_entries_to_df(
    entries: list[dict],
    pair_id: str,
    tumor_sample: str,
    normal_sample: str,
    status: str,
    reason: str = "",
) -> pd.DataFrame:
    columns = [
        "peptide_id",
        "pair_id",
        "tumor_sample",
        "normal_sample",
        "mhc_class",
        "peptide",
        "peptide_length",
        "peptide_qc_status",
        "peptide_qc_reasons",
        "peptide_qc_flags",
        "parent_record_count",
        "protein_count",
        "taxon_count",
    ]
    rows = []
    for idx, entry in enumerate(entries, start=1):
        rows.append(
            {
                "peptide_id": f"microbial_core_{idx:07d}" if status == "core" else "",
                "pair_id": pair_id,
                "tumor_sample": tumor_sample,
                "normal_sample": normal_sample,
                "mhc_class": entry["mhc_class"],
                "peptide": entry["peptide"],
                "peptide_length": entry["peptide_length"],
                "peptide_qc_status": status,
                "peptide_qc_reasons": entry.get("peptide_qc_reasons", reason),
                "peptide_qc_flags": entry.get("peptide_qc_flags", ""),
                "parent_record_count": len(entry["parent_record_ids"]),
                "protein_count": len(entry["protein_accessions"]),
                "taxon_count": len(entry["ctaxa_ids"]),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def peptide_blacklist_status(entry: dict) -> str:
    statuses = {str(s) for s in entry.get("blacklist_statuses", set()) if str(s)}
    if statuses and statuses == {"exclusive_contaminant"}:
        return "exclusive_contaminant"
    flags = []
    if "exclusive_contaminant" in statuses:
        flags.append("blacklist_partial_contaminant_retained")
    if "mixed" in statuses:
        flags.append("blacklist_mixed_retained")
    if "unresolved" in statuses:
        flags.append("blacklist_unresolved_retained")
    entry["peptide_qc_flags"] = ";".join(flags)
    return "retained"


def build_parent_map_df(
    core_entries: list[dict],
    excluded_entries: list[dict],
    parent_positions: dict[tuple[str, str, str], list[int]],
) -> pd.DataFrame:
    columns = [
        "mhc_class",
        "peptide",
        "peptide_length",
        "peptide_qc_status",
        "parent_record_id",
        "positions_1based",
    ]
    rows = []
    for status, entries in (("core", core_entries), ("excluded", excluded_entries)):
        for entry in entries:
            for parent_id in sorted(entry["parent_record_ids"]):
                positions = parent_positions.get((entry["mhc_class"], entry["peptide"], parent_id), [])
                rows.append(
                    {
                        "mhc_class": entry["mhc_class"],
                        "peptide": entry["peptide"],
                        "peptide_length": entry["peptide_length"],
                        "peptide_qc_status": status,
                        "parent_record_id": parent_id,
                        "positions_1based": ",".join(str(pos) for pos in sorted(set(positions))),
                    }
                )
    return pd.DataFrame(rows, columns=columns)


def build_matched_normal_df(
    excluded_entries: list[dict],
    normal_peptide_map: dict[tuple[str, str], dict],
) -> pd.DataFrame:
    columns = [
        "mhc_class",
        "peptide",
        "peptide_length",
        "normal_parent_record_count",
        "normal_protein_accessions",
        "normal_ctaxa_ids",
        "normal_parent_record_ids",
    ]
    rows = []
    for entry in excluded_entries:
        normal_entry = normal_peptide_map.get((entry["mhc_class"], entry["peptide"]), {})
        rows.append(
            {
                "mhc_class": entry["mhc_class"],
                "peptide": entry["peptide"],
                "peptide_length": entry["peptide_length"],
                "normal_parent_record_count": len(normal_entry.get("parent_record_ids", set())),
                "normal_protein_accessions": ",".join(sorted(normal_entry.get("protein_accessions", set()))),
                "normal_ctaxa_ids": ",".join(sorted(normal_entry.get("ctaxa_ids", set()))),
                "normal_parent_record_ids": ",".join(sorted(normal_entry.get("parent_record_ids", set()))),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def write_fasta(core_df: pd.DataFrame, mhc_class: str, path: Path) -> None:
    subset = core_df[core_df["mhc_class"] == mhc_class].copy()
    lines = []
    for row in subset.itertuples(index=False):
        lines.append(f">{row.peptide_id}|{row.mhc_class}|len={row.peptide_length}")
        lines.append(str(row.peptide))
    path.write_text("\n".join(lines) + ("\n" if lines else ""))


def write_combined_fasta(core_df: pd.DataFrame, path: Path) -> None:
    lines = []
    for row in core_df.itertuples(index=False):
        lines.append(f">{row.peptide_id}|{row.mhc_class}|len={row.peptide_length}")
        lines.append(str(row.peptide))
    path.write_text("\n".join(lines) + ("\n" if lines else ""))


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def write_stagewise_qc(path: Path, rows: list[tuple[str, int]]) -> None:
    pd.DataFrame([{"step": step, "n": int(count)} for step, count in rows]).to_csv(
        path,
        sep="\t",
        index=False,
    )


def build_pair_core(
    tumor_protein_hits: Path,
    normal_protein_hits: Path,
    outdir: Path,
    tumor_sample: str,
    normal_sample: str,
    pair_id: Optional[str] = None,
    blacklist: Optional[Path] = None,
    mhc_i_lengths: Optional[list[int]] = None,
    mhc_ii_lengths: Optional[list[int]] = None,
    min_pident: float = 100.0,
    max_evalue: float = 1e-5,
    min_qcovs: float = 90.0,
    allow_missing_blacklist: bool = False,
    blacklist_sha256: str = "",
) -> dict:
    mhc_i_lengths = mhc_i_lengths or [8, 9, 10, 11]
    mhc_ii_lengths = mhc_ii_lengths or [13, 14, 15, 16, 17]
    pair_id = pair_id or f"{tumor_sample},{normal_sample}"
    validate_pair_id(pair_id, tumor_sample, normal_sample)
    outdir.mkdir(parents=True, exist_ok=True)

    contaminant_taxids, verdict_by_taxid, actual_blacklist_sha256 = load_contaminant_taxids(
        blacklist,
        allow_missing=allow_missing_blacklist,
        expected_sha256=blacklist_sha256,
    )

    tumor_pass, tumor_excluded = read_protein_hits(
        tumor_protein_hits,
        tumor_sample,
        "tumor",
        min_pident,
        max_evalue,
        min_qcovs,
        contaminant_taxids,
        exclude_blacklist=False,
    )
    normal_pass, normal_excluded = read_protein_hits(
        normal_protein_hits,
        normal_sample,
        "normal",
        min_pident,
        max_evalue,
        min_qcovs,
        contaminant_taxids,
        exclude_blacklist=False,
    )

    normal_parent_hashes = set(normal_pass["parent_sequence_sha256"].dropna().astype(str))
    exact_parent_mask = tumor_pass["parent_sequence_sha256"].astype(str).isin(normal_parent_hashes)
    tumor_parent_exact_excluded = tumor_pass[exact_parent_mask].copy()
    if not tumor_parent_exact_excluded.empty:
        tumor_parent_exact_excluded["parent_qc_status"] = "excluded"
        tumor_parent_exact_excluded["parent_qc_reasons"] = "excluded_exact_parent_in_matched_normal"
    tumor_core = tumor_pass[~exact_parent_mask].copy().reset_index(drop=True)

    parent_excluded = pd.concat(
        [tumor_excluded, normal_excluded, tumor_parent_exact_excluded],
        ignore_index=True,
    )

    tumor_peptide_map, tumor_parent_positions = build_peptide_maps(
        tumor_core,
        mhc_i_lengths,
        mhc_ii_lengths,
    )
    normal_peptide_map, _ = build_peptide_maps(
        normal_pass,
        mhc_i_lengths,
        mhc_ii_lengths,
    )
    normal_peptide_keys = set(normal_peptide_map)

    core_entries = []
    matched_normal_excluded_entries = []
    blacklist_excluded_entries = []
    for key in sorted(tumor_peptide_map, key=lambda x: (x[0], len(x[1]), x[1])):
        entry = tumor_peptide_map[key]
        if key in normal_peptide_keys:
            entry["peptide_qc_reasons"] = "excluded_exact_peptide_in_matched_normal"
            matched_normal_excluded_entries.append(entry)
        elif peptide_blacklist_status(entry) == "exclusive_contaminant":
            entry["peptide_qc_reasons"] = "excluded_exclusive_contaminant"
            blacklist_excluded_entries.append(entry)
        else:
            core_entries.append(entry)
    excluded_entries = matched_normal_excluded_entries + blacklist_excluded_entries

    core_df = peptide_entries_to_df(
        core_entries,
        pair_id,
        tumor_sample,
        normal_sample,
        "core",
    )
    excluded_peptide_df = peptide_entries_to_df(
        excluded_entries,
        pair_id,
        tumor_sample,
        normal_sample,
        "excluded",
    )
    parent_map_df = build_parent_map_df(core_entries, excluded_entries, tumor_parent_positions)
    matched_normal_df = build_matched_normal_df(matched_normal_excluded_entries, normal_peptide_map)

    write_tsv(tumor_core, outdir / "microbial_parent_core.tsv")
    write_tsv(parent_excluded, outdir / "microbial_parent_excluded.tsv")
    write_tsv(core_df, outdir / "microbial_peptide_core.tsv")
    write_tsv(excluded_peptide_df, outdir / "microbial_peptide_excluded.tsv")
    write_tsv(parent_map_df, outdir / "microbial_peptide_parent_map.tsv")
    write_tsv(matched_normal_df, outdir / "matched_normal_peptide.tsv")
    write_combined_fasta(core_df, outdir / "microbial_peptide_core.fasta")
    write_fasta(core_df, MHC_CLASS_I, outdir / "microbial_peptide_core_hla_i.fasta")
    write_fasta(core_df, MHC_CLASS_II, outdir / "microbial_peptide_core_hla_ii.fasta")

    stagewise_rows = [
        ("tumor_parent_input_rows", len(tumor_pass) + len(tumor_excluded)),
        ("normal_parent_input_rows", len(normal_pass) + len(normal_excluded)),
        ("tumor_parent_qc_pass", len(tumor_pass)),
        ("normal_parent_qc_pass", len(normal_pass)),
        ("tumor_parent_qc_excluded", len(tumor_excluded)),
        ("normal_parent_qc_excluded", len(normal_excluded)),
        ("tumor_parent_exact_normal_excluded", len(tumor_parent_exact_excluded)),
        ("tumor_parent_core", len(tumor_core)),
        ("tumor_unique_peptides_before_normal_subtraction", len(tumor_peptide_map)),
        ("normal_unique_peptides_for_subtraction", len(normal_peptide_map)),
        ("tumor_unique_peptides_excluded_exact_normal", len(matched_normal_excluded_entries)),
        ("tumor_unique_peptides_excluded_blacklist", len(blacklist_excluded_entries)),
        ("tumor_unique_peptides_core", len(core_entries)),
        ("tumor_unique_peptides_core_hla_i", int((core_df["mhc_class"] == MHC_CLASS_I).sum()) if not core_df.empty else 0),
        ("tumor_unique_peptides_core_hla_ii", int((core_df["mhc_class"] == MHC_CLASS_II).sum()) if not core_df.empty else 0),
    ]
    write_stagewise_qc(outdir / "stagewise_qc.tsv", stagewise_rows)

    manifest = {
        "policy_version": POLICY_VERSION,
        "pair_id": pair_id,
        "tumor_sample": tumor_sample,
        "normal_sample": normal_sample,
        "tumor_protein_hits": str(tumor_protein_hits),
        "normal_protein_hits": str(normal_protein_hits),
        "tumor_protein_hits_sha256": sha256_file(tumor_protein_hits),
        "normal_protein_hits_sha256": sha256_file(normal_protein_hits),
        "blacklist": str(blacklist) if blacklist else "",
        "blacklist_sha256": actual_blacklist_sha256,
        "blacklist_sha256_expected": blacklist_sha256,
        "allow_missing_blacklist": allow_missing_blacklist,
        "blacklist_evaluation": "enabled" if blacklist else "not_evaluated",
        "blacklist_taxids": len(verdict_by_taxid),
        "contaminant_taxids": len(contaminant_taxids),
        "mhc_i_lengths": mhc_i_lengths,
        "mhc_ii_lengths": mhc_ii_lengths,
        "min_pident": min_pident,
        "max_evalue": max_evalue,
        "min_qcovs": min_qcovs,
        "outputs": {
            "microbial_parent_core": str(outdir / "microbial_parent_core.tsv"),
            "microbial_parent_excluded": str(outdir / "microbial_parent_excluded.tsv"),
            "microbial_peptide_core": str(outdir / "microbial_peptide_core.tsv"),
            "microbial_peptide_core_fasta": str(outdir / "microbial_peptide_core.fasta"),
            "microbial_peptide_core_hla_i_fasta": str(outdir / "microbial_peptide_core_hla_i.fasta"),
            "microbial_peptide_core_hla_ii_fasta": str(outdir / "microbial_peptide_core_hla_ii.fasta"),
            "microbial_peptide_excluded": str(outdir / "microbial_peptide_excluded.tsv"),
            "microbial_peptide_parent_map": str(outdir / "microbial_peptide_parent_map.tsv"),
            "matched_normal_peptide": str(outdir / "matched_normal_peptide.tsv"),
            "stagewise_qc": str(outdir / "stagewise_qc.tsv"),
        },
        "stagewise_counts": {step: int(count) for step, count in stagewise_rows},
    }
    with (outdir / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-id", default="")
    parser.add_argument("--tumor-sample", required=True)
    parser.add_argument("--normal-sample", required=True)
    parser.add_argument("--tumor-protein-hits", required=True)
    parser.add_argument("--normal-protein-hits", required=True)
    parser.add_argument("--blacklist", default="")
    parser.add_argument("--blacklist-sha256", default="")
    parser.add_argument(
        "--allow-missing-blacklist",
        action="store_true",
        help="Exploratory mode only: continue when no contaminant blacklist is configured.",
    )
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--mhc-i-lengths", default="8,9,10,11")
    parser.add_argument("--mhc-ii-lengths", default="13,14,15,16,17")
    parser.add_argument("--min-pident", type=float, default=100.0)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--min-qcovs", type=float, default=90.0)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = build_pair_core(
        tumor_protein_hits=Path(args.tumor_protein_hits),
        normal_protein_hits=Path(args.normal_protein_hits),
        outdir=Path(args.outdir),
        tumor_sample=args.tumor_sample,
        normal_sample=args.normal_sample,
        pair_id=args.pair_id or None,
        blacklist=Path(args.blacklist) if args.blacklist else None,
        mhc_i_lengths=parse_int_list(args.mhc_i_lengths),
        mhc_ii_lengths=parse_int_list(args.mhc_ii_lengths),
        min_pident=args.min_pident,
        max_evalue=args.max_evalue,
        min_qcovs=args.min_qcovs,
        allow_missing_blacklist=args.allow_missing_blacklist,
        blacklist_sha256=args.blacklist_sha256,
    )
    print(json.dumps(manifest["stagewise_counts"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
