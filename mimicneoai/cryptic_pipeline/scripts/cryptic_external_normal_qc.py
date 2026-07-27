#!/usr/bin/env python3
"""Apply frozen external-normal resources to a Cryptic peptide Core.

This stage consumes 08b-cryptic_core_qc peptide Core outputs and materializes
the tumor-restricted primary Core used by downstream binding. It uses exact
peptide sequence and HLA-class matching only.
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


POLICY_VERSION = "cryptic_external_normal_qc_v1.0"
EXPECTED_UPSTREAM_POLICY_VERSION = "cryptic_core_qc_v1.0"
EXPECTED_RESOURCE_FREEZE_VERSION = "external_normal_resources_v1.0_20260727"
COORDINATE_STATUS = "not_evaluable_coordinates_not_materialized"
SMORF_DETECTED = "normal_translation_detected_by_smorf_atlas"
SMORF_NOT_DETECTED = "normal_translation_not_detected_in_evaluated_smorf_atlas"
SMORF_NOT_EVALUATED = "normal_translation_not_evaluated"
HLA_DETECTED = "normal_hla_presentation_detected"
HLA_NOT_DETECTED = "normal_hla_presentation_not_detected_in_evaluated_hla_ligand_atlas"
HLA_NOT_EVALUATED = "normal_hla_presentation_not_evaluated"
RETAINED_STATUS = "not_detected_in_evaluated_frozen_resources"
EXCLUDED_STATUS = "external_normal_exact_match_excluded"
NOT_EVALUATED_STATUS = "external_normal_resources_not_evaluated"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
ALLOWED_LENGTHS = {
    "HLA-I": {8, 9, 10, 11},
    "HLA-II": {13, 14, 15, 16, 17},
}


def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


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


def read_json(path: Path) -> dict[str, object]:
    with path.open() as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def write_tsv(rows: list[dict[str, object]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(f">{row['peptide_record_id']}\n{row['peptide']}\n")


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


def normalize_hla_class(value: object) -> str:
    text = str(value).strip().upper().replace("_", "-")
    aliases = {
        "MHC-I": "HLA-I",
        "MHCI": "HLA-I",
        "CLASS-I": "HLA-I",
        "HLA-I": "HLA-I",
        "MHC-II": "HLA-II",
        "MHCII": "HLA-II",
        "CLASS-II": "HLA-II",
        "HLA-II": "HLA-II",
    }
    return aliases.get(text, text)


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)


def manifest_processed_entry(resource_manifest: dict[str, object], rel_path: str) -> dict[str, object]:
    processed = resource_manifest.get("processed_files", {})
    if not isinstance(processed, dict) or rel_path not in processed:
        raise ValueError(f"resource manifest missing processed_files entry: {rel_path}")
    entry = processed[rel_path]
    if not isinstance(entry, dict):
        raise ValueError(f"resource manifest entry is invalid: {rel_path}")
    return entry


def resolve_resource_path(manifest_path: Path, explicit: str, rel_path: str) -> Path:
    if explicit:
        return Path(explicit)
    return manifest_path.parent / rel_path


def validate_resource_file(path: Path, manifest_entry: dict[str, object]) -> dict[str, object]:
    identity = file_identity(path)
    if not identity.get("exists"):
        raise FileNotFoundError(f"external-normal resource file not found: {path}")
    expected_size = int(manifest_entry.get("size_bytes", -1))
    expected_sha = str(manifest_entry.get("sha256", ""))
    if identity.get("size") != expected_size or identity.get("sha256") != expected_sha:
        raise ValueError(
            f"external-normal resource identity mismatch for {path}: "
            f"expected size={expected_size} sha256={expected_sha}, "
            f"observed size={identity.get('size')} sha256={identity.get('sha256')}"
        )
    return identity


def validate_resource_contract(resource_manifest: dict[str, object]) -> None:
    freeze_version = str(resource_manifest.get("freeze_version", ""))
    if freeze_version != EXPECTED_RESOURCE_FREEZE_VERSION:
        raise ValueError(
            f"unsupported external-normal resource freeze_version: {freeze_version}; "
            f"expected {EXPECTED_RESOURCE_FREEZE_VERSION}"
        )
    policy = resource_manifest.get("matching_policy", {})
    if not isinstance(policy, dict):
        raise ValueError("resource manifest missing matching_policy")
    hla_i_lengths = [int(v) for v in policy.get("hla_i_lengths", [])]
    hla_ii_lengths = [int(v) for v in policy.get("hla_ii_lengths", [])]
    if hla_i_lengths != sorted(ALLOWED_LENGTHS["HLA-I"]):
        raise ValueError("resource manifest HLA-I length policy must be 8-11")
    if hla_ii_lengths != sorted(ALLOWED_LENGTHS["HLA-II"]):
        raise ValueError("resource manifest HLA-II length policy must be 13-17")
    smorf_policy = str(policy.get("smorf", "")).lower()
    hla_policy = str(policy.get("hla_ligand", "")).lower()
    if "exact peptide" not in smorf_policy:
        raise ValueError("resource manifest smORF policy must be exact peptide matching")
    if "exact peptide" not in hla_policy or "same hla class" not in hla_policy:
        raise ValueError("resource manifest HLA ligand policy must be exact peptide plus same HLA class")


def resource_identities(args: argparse.Namespace) -> tuple[bool, dict[str, object], dict[str, Path]]:
    manifest_path = Path(args.resource_manifest) if args.resource_manifest else None
    if manifest_path is None or not manifest_path.exists():
        if args.allow_missing_external_normal_resources:
            return False, {
                "resource_manifest": {
                    "path": str(manifest_path) if manifest_path is not None else "",
                    "exists": False,
                }
            }, {}
        raise FileNotFoundError("resource_manifest is required for formal CrypticExternalNormalQC")

    resource_manifest = read_json(manifest_path)
    validate_resource_contract(resource_manifest)
    rel_paths = {
        "smorf_match_index": "processed/normal_smorf_peptide_match_index.tsv.gz",
        "smorf_parent_map": "processed/normal_smorf_peptide_parent_map.tsv.gz",
        "hla_ligand_match_index": "processed/normal_hla_peptide_match_index.tsv.gz",
        "hla_ligand_evidence": "processed/normal_hla_ligand_evidence.tsv.gz",
    }
    explicit = {
        "smorf_match_index": args.smorf_match_index,
        "smorf_parent_map": args.smorf_parent_map,
        "hla_ligand_match_index": args.hla_ligand_match_index,
        "hla_ligand_evidence": args.hla_ligand_evidence,
    }
    identities: dict[str, object] = {
        "resource_manifest": file_identity(manifest_path),
        "resource_manifest_sha256": sha256_file(manifest_path),
        "freeze_version": resource_manifest.get("freeze_version", ""),
    }
    paths: dict[str, Path] = {}
    try:
        for key, rel_path in rel_paths.items():
            path = resolve_resource_path(manifest_path, str(explicit[key] or ""), rel_path)
            paths[key] = path
            identities[key] = validate_resource_file(path, manifest_processed_entry(resource_manifest, rel_path))
    except Exception:
        if args.allow_missing_external_normal_resources:
            return False, identities, paths
        raise
    return True, identities, paths


def validate_signature_entry(actual_path: Path, entry: object, label: str) -> dict[str, object]:
    if not isinstance(entry, dict):
        raise ValueError(f"upstream output_signature entry for {label} is invalid")
    actual = file_identity(actual_path)
    if not actual.get("exists"):
        raise FileNotFoundError(f"required upstream 08b output does not exist: {actual_path}")
    if not bool(entry.get("exists", False)):
        raise ValueError(f"upstream output_signature says {label} does not exist")
    if actual.get("size") != entry.get("size") or actual.get("sha256") != entry.get("sha256"):
        raise ValueError(
            f"upstream output_signature mismatch for {label}: "
            f"expected size={entry.get('size')} sha256={entry.get('sha256')}, "
            f"observed size={actual.get('size')} sha256={actual.get('sha256')}"
        )
    return actual


def validate_upstream_08b(core_path: Path, parent_map_path: Path, manifest_path: Path) -> dict[str, object]:
    for label, path in [
        ("cryptic_peptide_core.tsv", core_path),
        ("cryptic_peptide_parent_map.tsv", parent_map_path),
        ("upstream run_manifest.json", manifest_path),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"required upstream 08b input missing: {label}: {path}")
    manifest = read_json(manifest_path)
    if manifest.get("run_status") != "complete":
        raise ValueError("upstream 08b run_manifest.json must have run_status=complete")
    if manifest.get("policy_version") != EXPECTED_UPSTREAM_POLICY_VERSION:
        raise ValueError(
            f"upstream 08b policy_version must be {EXPECTED_UPSTREAM_POLICY_VERSION}; "
            f"observed {manifest.get('policy_version', '')}"
        )
    output_signature = manifest.get("output_signature", {})
    if not isinstance(output_signature, dict):
        raise ValueError("upstream 08b run_manifest.json missing output_signature")
    validate_signature_entry(core_path, output_signature.get("cryptic_peptide_core.tsv"), "cryptic_peptide_core.tsv")
    validate_signature_entry(
        parent_map_path,
        output_signature.get("cryptic_peptide_parent_map.tsv"),
        "cryptic_peptide_parent_map.tsv",
    )
    return manifest


def add_key_columns(df: pd.DataFrame, peptide_col: str = "peptide_sequence") -> pd.DataFrame:
    out = df.copy()
    out["__hla_class"] = out["hla_class"].map(normalize_hla_class)
    out["__peptide"] = out[peptide_col].astype(str).str.upper()
    return out


def build_lookup(df: pd.DataFrame, peptide_col: str = "peptide_sequence") -> set[tuple[str, str]]:
    keyed = add_key_columns(df, peptide_col=peptide_col)
    return set(zip(keyed["__hla_class"], keyed["__peptide"]))


def group_rows_by_key(df: pd.DataFrame, peptide_col: str = "peptide_sequence") -> dict[tuple[str, str], list[dict[str, str]]]:
    keyed = add_key_columns(df, peptide_col=peptide_col)
    grouped: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in keyed.to_dict("records"):
        grouped[(str(row["__hla_class"]), str(row["__peptide"]))].append(row)
    return dict(grouped)


def required_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} missing required columns: {', '.join(missing)}")


def parse_int(value: object, label: str) -> int:
    try:
        return int(str(value))
    except Exception as exc:
        raise ValueError(f"{label} must be an integer: {value}") from exc


def validate_peptide_row(row: dict[str, object], sample: str, table_label: str, row_id: str) -> None:
    observed_sample = str(row.get("sample", ""))
    if observed_sample != sample:
        raise ValueError(f"{table_label} sample mismatch for {row_id}: {observed_sample} != {sample}")
    hla_class = normalize_hla_class(row.get("mhc_class", ""))
    if hla_class not in ALLOWED_LENGTHS:
        raise ValueError(f"{table_label} has unsupported mhc_class for {row_id}: {row.get('mhc_class', '')}")
    peptide = str(row.get("peptide", ""))
    if not peptide or peptide.upper() != peptide or any(residue not in STANDARD_AA for residue in peptide):
        raise ValueError(f"{table_label} has non-standard peptide sequence for {row_id}: {peptide}")
    peptide_length = parse_int(row.get("peptide_length", ""), f"{table_label} peptide_length for {row_id}")
    if peptide_length != len(peptide):
        raise ValueError(
            f"{table_label} peptide_length does not match sequence length for {row_id}: "
            f"{peptide_length} != {len(peptide)}"
        )
    if peptide_length not in ALLOWED_LENGTHS[hla_class]:
        allowed = "-".join([str(min(ALLOWED_LENGTHS[hla_class])), str(max(ALLOWED_LENGTHS[hla_class]))])
        raise ValueError(f"{table_label} {hla_class} peptide length for {row_id} must be {allowed}: {peptide_length}")


def validate_core_parent_contract(core: pd.DataFrame, parent_map: pd.DataFrame, sample: str) -> None:
    if core["peptide_record_id"].duplicated().any():
        duplicates = core.loc[core["peptide_record_id"].duplicated(), "peptide_record_id"].head(5).tolist()
        raise ValueError(f"cryptic_peptide_core.tsv must contain unique peptide_record_id values: {duplicates}")

    core_records = {str(row["peptide_record_id"]): row for row in core.to_dict("records")}
    parent_records = parent_map.to_dict("records")
    parent_ids = {str(row.get("peptide_record_id", "")) for row in parent_records}
    core_ids = set(core_records)

    orphan_ids = sorted(parent_ids - core_ids)
    if orphan_ids:
        raise ValueError(f"cryptic_peptide_parent_map.tsv has orphan peptide_record_id values: {orphan_ids[:5]}")
    missing_parent_ids = sorted(core_ids - parent_ids)
    if missing_parent_ids:
        raise ValueError(f"cryptic_peptide_core.tsv peptides missing parent-map records: {missing_parent_ids[:5]}")

    for peptide_id, row in core_records.items():
        validate_peptide_row(row, sample, "cryptic_peptide_core.tsv", peptide_id)

    for idx, row in enumerate(parent_records, start=1):
        peptide_id = str(row.get("peptide_record_id", ""))
        validate_peptide_row(row, sample, "cryptic_peptide_parent_map.tsv", peptide_id or f"row{idx}")
        core_row = core_records.get(peptide_id)
        if core_row is None:
            continue
        mismatches = []
        for column in ["sample", "peptide", "peptide_length"]:
            if str(row.get(column, "")) != str(core_row.get(column, "")):
                mismatches.append(column)
        if normalize_hla_class(row.get("mhc_class", "")) != normalize_hla_class(core_row.get("mhc_class", "")):
            mismatches.append("mhc_class")
        if mismatches:
            raise ValueError(
                f"Core and parent-map disagree for peptide_record_id {peptide_id}: {', '.join(mismatches)}"
            )


def select_core_parent_map_rows(parent_map: pd.DataFrame, core_ids: set[str]) -> pd.DataFrame:
    if "peptide_core_status" not in parent_map.columns:
        return parent_map
    status = parent_map["peptide_core_status"].astype(str)
    core_parent_map = parent_map[status.eq("core")].copy()
    noncore_ids = {
        str(value)
        for value in parent_map.loc[~status.eq("core"), "peptide_record_id"].tolist()
        if str(value).strip()
    }
    conflicting_ids = sorted(noncore_ids & core_ids)
    if conflicting_ids:
        raise ValueError(f"non-core parent-map rows reference Core peptide_record_id values: {conflicting_ids[:5]}")
    return core_parent_map


def compact_values(rows: list[dict[str, str]], key: str) -> str:
    vals = sorted({str(row.get(key, "")).strip() for row in rows if str(row.get(key, "")).strip()})
    return ";".join(vals)


def base_output_row(row: dict[str, object], resource_freeze_version: str) -> dict[str, object]:
    out = dict(row)
    out["hla_class_normalized"] = normalize_hla_class(row.get("mhc_class", ""))
    out["normal_smorf_coordinate_status"] = COORDINATE_STATUS
    out["resource_freeze_version"] = resource_freeze_version
    return out


def make_source_evidence(row: dict[str, object], resource_freeze_version: str) -> dict[str, object]:
    out = base_output_row(row, resource_freeze_version)
    out.update({
        "evidence_type": "source_core_parent_map",
        "source_resource": "08b-cryptic_core_qc",
        "resource_record_id": str(row.get("parent_record_id", "")),
        "normal_translation_status": str(row.get("normal_translation_status", "")),
        "normal_hla_presentation_status": str(row.get("normal_hla_presentation_status", "")),
        "external_normal_status": str(row.get("external_normal_status", "")),
        "smorf_iorf_id": "",
        "smorf_orf_id": "",
        "smorf_peptide_start_1based": "",
        "hla_ligand_peptide_sequence_id": "",
        "hla_ligand_donor": "",
        "hla_ligand_tissue": "",
        "hla_ligand_donor_hla_genotype": "",
    })
    return out


def make_smorf_evidence(row: dict[str, object], evidence: dict[str, str], resource_freeze_version: str) -> dict[str, object]:
    out = base_output_row(row, resource_freeze_version)
    source_record_id = str(evidence.get("source_record_id", ""))
    out.update({
        "evidence_type": "normal_smorf_translation",
        "source_resource": str(evidence.get("source_resource", "Chothani_2022_smORF_atlas")),
        "resource_record_id": source_record_id or str(evidence.get("iorf_id", "")),
        "normal_translation_status": SMORF_DETECTED,
        "normal_hla_presentation_status": str(row.get("normal_hla_presentation_status", "")),
        "external_normal_status": EXCLUDED_STATUS,
        "smorf_iorf_id": str(evidence.get("iorf_id", "")),
        "smorf_orf_id": str(evidence.get("orf_id", "")),
        "smorf_peptide_start_1based": str(evidence.get("peptide_start_1based", "")),
        "hla_ligand_peptide_sequence_id": "",
        "hla_ligand_donor": "",
        "hla_ligand_tissue": "",
        "hla_ligand_donor_hla_genotype": "",
    })
    return out


def make_hla_evidence(row: dict[str, object], evidence: dict[str, str], resource_freeze_version: str) -> dict[str, object]:
    out = base_output_row(row, resource_freeze_version)
    donor = str(evidence.get("donor", ""))
    tissue = str(evidence.get("tissue", ""))
    seq_id = str(evidence.get("peptide_sequence_id", ""))
    out.update({
        "evidence_type": "normal_hla_ligand",
        "source_resource": str(evidence.get("source_resource", "HLA_Ligand_Atlas_2020.12")),
        "resource_record_id": "|".join(token for token in [seq_id, donor, tissue] if token),
        "normal_translation_status": str(row.get("normal_translation_status", "")),
        "normal_hla_presentation_status": HLA_DETECTED,
        "external_normal_status": EXCLUDED_STATUS,
        "smorf_iorf_id": "",
        "smorf_orf_id": "",
        "smorf_peptide_start_1based": "",
        "hla_ligand_peptide_sequence_id": seq_id,
        "hla_ligand_donor": donor,
        "hla_ligand_tissue": tissue,
        "hla_ligand_donor_hla_genotype": str(evidence.get("donor_hla_genotype", "")),
    })
    return out


def build_external_normal_qc(args: argparse.Namespace) -> dict[str, object]:
    started_at = ts()
    if args.coordinate_matching_enabled:
        raise ValueError("coordinate matching is not implemented in policy v1.0")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    script_sha256 = sha256_file(Path(__file__))

    core_path = Path(args.cryptic_peptide_core)
    parent_map_path = Path(args.cryptic_peptide_parent_map)
    upstream_manifest_path = Path(args.upstream_manifest)

    upstream_manifest = validate_upstream_08b(core_path, parent_map_path, upstream_manifest_path)
    resources_evaluated, resource_identity, resource_paths = resource_identities(args)
    freeze_version = str(resource_identity.get("freeze_version", "")) if resources_evaluated else ""

    input_signature = {
        "policy_version": POLICY_VERSION,
        "script_sha256": script_sha256,
        "sample": args.sample,
        "allow_missing_external_normal_resources": bool(args.allow_missing_external_normal_resources),
        "coordinate_matching_enabled": bool(args.coordinate_matching_enabled),
        "hla_class_normalization": "MHC-I/HLA-I->HLA-I; MHC-II/HLA-II->HLA-II",
        "exact_match_policy": "exact peptide sequence plus normalized HLA class",
        "coordinate_matching_status": COORDINATE_STATUS,
        "upstream_inputs": {
            "cryptic_peptide_core": file_identity(core_path),
            "cryptic_peptide_parent_map": file_identity(parent_map_path),
            "upstream_manifest": file_identity(upstream_manifest_path),
        },
        "upstream_08b": {
            "policy_version": upstream_manifest.get("policy_version", ""),
            "run_status": upstream_manifest.get("run_status", ""),
            "cryptic_peptide_core_output_signature": upstream_manifest.get("output_signature", {}).get(
                "cryptic_peptide_core.tsv", {}
            ),
            "cryptic_peptide_parent_map_output_signature": upstream_manifest.get("output_signature", {}).get(
                "cryptic_peptide_parent_map.tsv", {}
            ),
        },
        "resource_identities": resource_identity,
    }

    output_paths = [
        outdir / "cryptic_external_normal_annotated_core.tsv",
        outdir / "cryptic_external_normal_evidence.tsv",
        outdir / "cryptic_external_normal_excluded.tsv",
        outdir / "cryptic_tumor_restricted_primary_core.tsv",
        outdir / "cryptic_tumor_restricted_primary_core_hla_i.tsv",
        outdir / "cryptic_tumor_restricted_primary_core_hla_ii.tsv",
        outdir / "cryptic_tumor_restricted_primary_core_hla_i.fasta",
        outdir / "cryptic_tumor_restricted_primary_core_hla_ii.fasta",
        outdir / "cryptic_tumor_restricted_primary_core.fasta",
        outdir / "stagewise_qc.tsv",
    ]
    manifest_path = outdir / "run_manifest.json"
    if manifest_path.exists():
        old = read_json(manifest_path)
        if old.get("input_signature") != input_signature:
            raise ValueError(
                f"Existing CrypticExternalNormalQC manifest does not match current inputs/config: {manifest_path}"
            )
        old_output_signature = old.get("output_signature")
        if isinstance(old_output_signature, dict) and outputs_match_signature(old_output_signature):
            return {**old, "reused": True}

    core = load_table(core_path)
    parent_map = load_table(parent_map_path)
    required_columns(
        core,
        ["sample", "peptide_record_id", "mhc_class", "peptide_length", "peptide"],
        "cryptic_peptide_core.tsv",
    )
    required_columns(
        parent_map,
        ["sample", "peptide_record_id", "parent_record_id", "source_parent_id", "mhc_class", "peptide_length", "peptide"],
        "cryptic_peptide_parent_map.tsv",
    )
    if core["peptide_record_id"].duplicated().any():
        raise ValueError("cryptic_peptide_core.tsv must contain unique peptide_record_id values")
    core_ids = {str(value) for value in core["peptide_record_id"].tolist()}
    core_parent_map = select_core_parent_map_rows(parent_map, core_ids)
    validate_core_parent_contract(core, core_parent_map, args.sample)

    smorf_lookup: set[tuple[str, str]] = set()
    hla_lookup: set[tuple[str, str]] = set()
    smorf_evidence: dict[tuple[str, str], list[dict[str, str]]] = {}
    hla_evidence: dict[tuple[str, str], list[dict[str, str]]] = {}
    if resources_evaluated and not core.empty:
        smorf_index = load_table(resource_paths["smorf_match_index"])
        hla_index = load_table(resource_paths["hla_ligand_match_index"])
        smorf_parent_map = load_table(resource_paths["smorf_parent_map"])
        hla_ligand_evidence = load_table(resource_paths["hla_ligand_evidence"])
        required_columns(smorf_index, ["hla_class", "peptide_sequence"], "normal_smorf_peptide_match_index")
        required_columns(hla_index, ["hla_class", "peptide_sequence"], "normal_hla_peptide_match_index")
        required_columns(smorf_parent_map, ["hla_class", "peptide_sequence"], "normal_smorf_peptide_parent_map")
        required_columns(hla_ligand_evidence, ["hla_class", "peptide_sequence"], "normal_hla_ligand_evidence")
        smorf_lookup = build_lookup(smorf_index)
        hla_lookup = build_lookup(hla_index)
        smorf_evidence = group_rows_by_key(smorf_parent_map)
        hla_evidence = group_rows_by_key(hla_ligand_evidence)

    core_records: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []
    retained_rows: list[dict[str, object]] = []
    smorf_match_count = 0
    hla_match_count = 0
    both_match_count = 0

    for record in core.to_dict("records"):
        row = base_output_row(record, freeze_version)
        row["normal_smorf_coordinate_status"] = COORDINATE_STATUS
        key = (str(row["hla_class_normalized"]), str(row["peptide"]).upper())
        smorf_match = resources_evaluated and key in smorf_lookup
        hla_match = resources_evaluated and key in hla_lookup

        if resources_evaluated:
            row["normal_translation_status"] = SMORF_DETECTED if smorf_match else SMORF_NOT_DETECTED
            row["normal_hla_presentation_status"] = HLA_DETECTED if hla_match else HLA_NOT_DETECTED
            row["normal_smorf_support_count"] = len(smorf_evidence.get(key, []))
            row["normal_hla_ligand_support_count"] = len(hla_evidence.get(key, []))
            if smorf_match or hla_match:
                reasons = []
                if smorf_match:
                    reasons.append("normal_smorf_exact_peptide_match")
                    smorf_match_count += 1
                if hla_match:
                    reasons.append("normal_hla_ligand_exact_peptide_match")
                    hla_match_count += 1
                if smorf_match and hla_match:
                    both_match_count += 1
                row["external_normal_status"] = EXCLUDED_STATUS
                row["external_normal_qc_reasons"] = ";".join(reasons)
                excluded_rows.append(row)
            else:
                row["external_normal_status"] = RETAINED_STATUS
                row["external_normal_qc_reasons"] = ""
                retained_rows.append(row)
        else:
            row["normal_translation_status"] = SMORF_NOT_EVALUATED
            row["normal_hla_presentation_status"] = HLA_NOT_EVALUATED
            row["normal_smorf_support_count"] = 0
            row["normal_hla_ligand_support_count"] = 0
            row["external_normal_status"] = NOT_EVALUATED_STATUS
            row["external_normal_qc_reasons"] = "external_normal_resources_not_evaluated"

        row["normal_smorf_iorf_ids"] = compact_values(smorf_evidence.get(key, []), "iorf_id") if smorf_match else ""
        row["normal_smorf_orf_ids"] = compact_values(smorf_evidence.get(key, []), "orf_id") if smorf_match else ""
        row["normal_hla_ligand_donors"] = compact_values(hla_evidence.get(key, []), "donor") if hla_match else ""
        row["normal_hla_ligand_tissues"] = compact_values(hla_evidence.get(key, []), "tissue") if hla_match else ""
        row["normal_hla_ligand_donor_hla_genotypes"] = (
            compact_values(hla_evidence.get(key, []), "donor_hla_genotype") if hla_match else ""
        )
        core_records.append(row)
        if smorf_match:
            for evidence in smorf_evidence.get(key, []):
                evidence_rows.append(make_smorf_evidence(row, evidence, freeze_version))
        if hla_match:
            for evidence in hla_evidence.get(key, []):
                evidence_rows.append(make_hla_evidence(row, evidence, freeze_version))

    core_by_peptide_id = {str(record.get("peptide_record_id", "")): record for record in core_records}
    for row in core_parent_map.to_dict("records"):
        peptide_id = str(row.get("peptide_record_id", ""))
        enriched = core_by_peptide_id.get(peptide_id, {})
        source_row = dict(row)
        for key in [
            "normal_translation_status",
            "normal_hla_presentation_status",
            "normal_smorf_coordinate_status",
            "external_normal_status",
        ]:
            source_row[key] = enriched.get(key, "")
        evidence_rows.append(make_source_evidence(source_row, freeze_version))

    if not resources_evaluated:
        retained_rows = []
        excluded_rows = list(core_records)

    hla_i_rows = [row for row in retained_rows if normalize_hla_class(row.get("mhc_class", "")) == "HLA-I"]
    hla_ii_rows = [row for row in retained_rows if normalize_hla_class(row.get("mhc_class", "")) == "HLA-II"]

    common_fields = list(core.columns)
    extra_fields = [
        "hla_class_normalized",
        "normal_translation_status",
        "normal_hla_presentation_status",
        "normal_smorf_coordinate_status",
        "normal_smorf_support_count",
        "normal_hla_ligand_support_count",
        "normal_smorf_iorf_ids",
        "normal_smorf_orf_ids",
        "normal_hla_ligand_donors",
        "normal_hla_ligand_tissues",
        "normal_hla_ligand_donor_hla_genotypes",
        "external_normal_status",
        "external_normal_qc_reasons",
        "resource_freeze_version",
    ]
    main_fields = common_fields + [field for field in extra_fields if field not in common_fields]
    evidence_fields = main_fields + [
        "evidence_type",
        "source_resource",
        "resource_record_id",
        "smorf_iorf_id",
        "smorf_orf_id",
        "smorf_peptide_start_1based",
        "hla_ligand_peptide_sequence_id",
        "hla_ligand_donor",
        "hla_ligand_tissue",
        "hla_ligand_donor_hla_genotype",
    ]

    if len(core_records) != len(core):
        raise AssertionError("annotated Core row count must equal input Core row count")
    if len(retained_rows) + len(excluded_rows) != len(core_records):
        raise AssertionError("retained + excluded rows must equal annotated Core rows")

    write_tsv(core_records, outdir / "cryptic_external_normal_annotated_core.tsv", main_fields)
    write_tsv(evidence_rows, outdir / "cryptic_external_normal_evidence.tsv", evidence_fields)
    write_tsv(excluded_rows, outdir / "cryptic_external_normal_excluded.tsv", main_fields)
    write_tsv(retained_rows, outdir / "cryptic_tumor_restricted_primary_core.tsv", main_fields)
    write_tsv(hla_i_rows, outdir / "cryptic_tumor_restricted_primary_core_hla_i.tsv", main_fields)
    write_tsv(hla_ii_rows, outdir / "cryptic_tumor_restricted_primary_core_hla_ii.tsv", main_fields)
    write_fasta(hla_i_rows, outdir / "cryptic_tumor_restricted_primary_core_hla_i.fasta")
    write_fasta(hla_ii_rows, outdir / "cryptic_tumor_restricted_primary_core_hla_ii.fasta")
    write_fasta(retained_rows, outdir / "cryptic_tumor_restricted_primary_core.fasta")

    stage_counts = [
        {"stage": "source_core_unique_peptides", "count": len(core_records)},
        {"stage": "annotated_core_unique_peptides", "count": len(core_records)},
        {"stage": "external_normal_resources_evaluated", "count": int(resources_evaluated)},
        {"stage": "normal_smorf_exact_match_peptides", "count": smorf_match_count},
        {"stage": "normal_hla_ligand_exact_match_peptides", "count": hla_match_count},
        {"stage": "both_external_normal_evidence_peptides", "count": both_match_count},
        {"stage": "external_normal_excluded_unique_peptides", "count": len(excluded_rows)},
        {
            "stage": "external_normal_not_evaluated_unique_peptides",
            "count": len(core_records) if not resources_evaluated else 0,
        },
        {"stage": "tumor_restricted_primary_core_unique_peptides", "count": len(retained_rows)},
        {"stage": "hla_i_tumor_restricted_primary_core_unique_peptides", "count": len(hla_i_rows)},
        {"stage": "hla_ii_tumor_restricted_primary_core_unique_peptides", "count": len(hla_ii_rows)},
        {"stage": "external_normal_evidence_rows", "count": len(evidence_rows)},
    ]
    write_tsv(stage_counts, outdir / "stagewise_qc.tsv", ["stage", "count"])

    outputs_hash = output_signature(output_paths)
    manifest = {
        "sample": args.sample,
        "policy_version": POLICY_VERSION,
        "run_status": "complete",
        "started_at": started_at,
        "finished_at": ts(),
        "code_commit": git_commit(),
        "script_sha256": script_sha256,
        "input_signature": input_signature,
        "upstream_08b_policy_version": upstream_manifest.get("policy_version", ""),
        "external_normal_evaluation": "evaluated" if resources_evaluated else "not_evaluated",
        "resource_freeze_version": freeze_version,
        "coordinate_matching_status": COORDINATE_STATUS,
        "outputs": {path.name: str(path) for path in output_paths},
        "output_signature": outputs_hash,
        "stage_counts": {row["stage"]: row["count"] for row in stage_counts},
        "binding_input_fasta": str(outdir / "cryptic_tumor_restricted_primary_core.fasta") if resources_evaluated else "",
        "binding_input_mode": "peptide-core" if resources_evaluated else "not_available",
        "reused": False,
    }
    write_json(manifest_path, manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", required=True)
    parser.add_argument("--cryptic-peptide-core", required=True)
    parser.add_argument("--cryptic-peptide-parent-map", required=True)
    parser.add_argument("--upstream-manifest", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--resource-manifest", default="")
    parser.add_argument("--smorf-match-index", default="")
    parser.add_argument("--smorf-parent-map", default="")
    parser.add_argument("--hla-ligand-match-index", default="")
    parser.add_argument("--hla-ligand-evidence", default="")
    parser.add_argument("--allow-missing-external-normal-resources", action="store_true")
    parser.add_argument("--coordinate-matching-enabled", action="store_true")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    manifest = build_external_normal_qc(build_parser().parse_args(argv))
    print(json.dumps(manifest, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
