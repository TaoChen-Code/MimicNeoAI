from __future__ import annotations

import argparse
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from mimicneoai.cryptic_pipeline.scripts.cryptic_external_normal_qc import (
    build_external_normal_qc,
    sha256_file,
)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def write_tsv_gz(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        pd.DataFrame(rows, columns=columns).to_csv(handle, sep="\t", index=False)


class CrypticExternalNormalQCTest(unittest.TestCase):
    core_columns = ["sample", "peptide_record_id", "mhc_class", "peptide_length", "peptide"]
    parent_columns = [
        "sample",
        "peptide_record_id",
        "parent_record_id",
        "source_parent_id",
        "mhc_class",
        "peptide_length",
        "peptide",
    ]

    def default_core_rows(self) -> list[dict[str, object]]:
        return [
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep1",
                "mhc_class": "MHC-I",
                "peptide_length": 8,
                "peptide": "ACDEFGHI",
            },
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep2",
                "mhc_class": "MHC-I",
                "peptide_length": 8,
                "peptide": "CDEFGHIK",
            },
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep3",
                "mhc_class": "MHC-II",
                "peptide_length": 13,
                "peptide": "ACDEFGHIKLMNP",
            },
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep4",
                "mhc_class": "MHC-I",
                "peptide_length": 8,
                "peptide": "KLMNPQRS",
            },
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep5",
                "mhc_class": "MHC-I",
                "peptide_length": 8,
                "peptide": "QRSTVWYA",
            },
        ]

    def default_parent_rows(self, core_rows: list[dict[str, object]]) -> list[dict[str, object]]:
        rows = [
            {
                "sample": row["sample"],
                "peptide_record_id": row["peptide_record_id"],
                "parent_record_id": f"parent_{row['peptide_record_id']}",
                "source_parent_id": f"source_{row['peptide_record_id']}",
                "mhc_class": row["mhc_class"],
                "peptide_length": row["peptide_length"],
                "peptide": row["peptide"],
            }
            for row in core_rows
        ]
        rows.append(
            {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "pep5",
                "parent_record_id": "parent_pep5_b",
                "source_parent_id": "source_pep5_b",
                "mhc_class": "MHC-I",
                "peptide_length": 8,
                "peptide": "QRSTVWYA",
            }
        )
        return rows

    def write_08b_manifest(
        self,
        manifest: Path,
        core: Path,
        parent_map: Path,
        *,
        run_status: str = "complete",
        policy_version: str = "cryptic_core_qc_v1.0",
        output_signature_override: dict[str, object] | None = None,
    ) -> None:
        output_signature = {
            "cryptic_peptide_core.tsv": {
                "path": str(core),
                "exists": core.exists(),
                "size": core.stat().st_size if core.exists() else 0,
                "sha256": sha256_file(core) if core.exists() else "",
            },
            "cryptic_peptide_parent_map.tsv": {
                "path": str(parent_map),
                "exists": parent_map.exists(),
                "size": parent_map.stat().st_size if parent_map.exists() else 0,
                "sha256": sha256_file(parent_map) if parent_map.exists() else "",
            },
        }
        if output_signature_override:
            output_signature.update(output_signature_override)
        manifest.write_text(
            json.dumps(
                {
                    "policy_version": policy_version,
                    "run_status": run_status,
                    "output_signature": output_signature,
                },
                indent=2,
            )
        )

    def _write_08b_core(
        self,
        root: Path,
        *,
        core_rows: list[dict[str, object]] | None = None,
        parent_rows: list[dict[str, object]] | None = None,
        run_status: str = "complete",
        policy_version: str = "cryptic_core_qc_v1.0",
    ) -> dict[str, Path]:
        core_rows = self.default_core_rows() if core_rows is None else core_rows
        parent_rows = self.default_parent_rows(core_rows) if parent_rows is None else parent_rows
        core = root / "08b" / "cryptic_peptide_core.tsv"
        parent_map = root / "08b" / "cryptic_peptide_parent_map.tsv"
        manifest = root / "08b" / "run_manifest.json"
        write_tsv(core, core_rows, self.core_columns)
        write_tsv(parent_map, parent_rows, self.parent_columns)
        self.write_08b_manifest(
            manifest,
            core,
            parent_map,
            run_status=run_status,
            policy_version=policy_version,
        )
        return {"core": core, "parent_map": parent_map, "manifest": manifest}

    def _write_resources(
        self,
        root: Path,
        *,
        freeze_version: str = "external_normal_resources_v1.0_20260727",
        hla_i_lengths: list[int] | None = None,
        hla_ii_lengths: list[int] | None = None,
    ) -> dict[str, Path]:
        resource_root = root / "resources"
        processed = resource_root / "processed"
        smorf_index = processed / "normal_smorf_peptide_match_index.tsv.gz"
        smorf_parent_map = processed / "normal_smorf_peptide_parent_map.tsv.gz"
        hla_index = processed / "normal_hla_peptide_match_index.tsv.gz"
        hla_evidence = processed / "normal_hla_ligand_evidence.tsv.gz"

        write_tsv_gz(
            smorf_index,
            [
                {"hla_class": "HLA-I", "peptide_sequence": "ACDEFGHI", "peptide_length": 8},
                {"hla_class": "HLA-II", "peptide_sequence": "CDEFGHIK", "peptide_length": 8},
                {"hla_class": "HLA-I", "peptide_sequence": "KLMNPQRS", "peptide_length": 8},
            ],
            ["hla_class", "peptide_sequence", "peptide_length"],
        )
        write_tsv_gz(
            smorf_parent_map,
            [
                {
                    "hla_class": "HLA-I",
                    "peptide_sequence": "ACDEFGHI",
                    "iorf_id": "iORF1",
                    "orf_id": "ORF1",
                    "source_record_id": "S1",
                    "peptide_start_1based": 1,
                },
                {
                    "hla_class": "HLA-II",
                    "peptide_sequence": "CDEFGHIK",
                    "iorf_id": "iORF_CLASS_MISMATCH",
                    "orf_id": "ORF_CLASS_MISMATCH",
                    "source_record_id": "S2",
                    "peptide_start_1based": 1,
                },
                {
                    "hla_class": "HLA-I",
                    "peptide_sequence": "KLMNPQRS",
                    "iorf_id": "iORF4",
                    "orf_id": "ORF4",
                    "source_record_id": "S4",
                    "peptide_start_1based": 2,
                },
            ],
            ["hla_class", "peptide_sequence", "iorf_id", "orf_id", "source_record_id", "peptide_start_1based"],
        )
        write_tsv_gz(
            hla_index,
            [
                {"hla_class": "HLA-II", "peptide_sequence": "ACDEFGHIKLMNP", "peptide_length": 13},
                {"hla_class": "HLA-I", "peptide_sequence": "KLMNPQRS", "peptide_length": 8},
            ],
            ["hla_class", "peptide_sequence", "peptide_length"],
        )
        write_tsv_gz(
            hla_evidence,
            [
                {
                    "hla_class": "HLA-II",
                    "peptide_sequence": "ACDEFGHIKLMNP",
                    "peptide_sequence_id": "HLA_SEQ_3",
                    "donor": "D1",
                    "donor_hla_genotype": "HLA-DRB1*01:01",
                    "tissue": "colon",
                    "source_resource": "HLA_Ligand_Atlas",
                },
                {
                    "hla_class": "HLA-I",
                    "peptide_sequence": "KLMNPQRS",
                    "peptide_sequence_id": "HLA_SEQ_4",
                    "donor": "D2",
                    "donor_hla_genotype": "HLA-A*02:01",
                    "tissue": "ileum",
                    "source_resource": "HLA_Ligand_Atlas",
                },
            ],
            [
                "hla_class",
                "peptide_sequence",
                "peptide_sequence_id",
                "donor",
                "donor_hla_genotype",
                "tissue",
                "source_resource",
            ],
        )

        manifest = resource_root / "resource_manifest.json"
        processed_files = {}
        for path in [smorf_index, smorf_parent_map, hla_index, hla_evidence]:
            rel = path.relative_to(resource_root).as_posix()
            processed_files[rel] = {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        manifest.write_text(
            json.dumps(
                {
                    "freeze_version": freeze_version,
                    "matching_policy": {
                        "smorf": "exact peptide sequence against windows generated from published smORFs",
                        "hla_ligand": "exact peptide sequence in the same HLA class",
                        "hla_i_lengths": hla_i_lengths or [8, 9, 10, 11],
                        "hla_ii_lengths": hla_ii_lengths or [13, 14, 15, 16, 17],
                    },
                    "processed_files": processed_files,
                },
                indent=2,
            )
        )
        return {
            "manifest": manifest,
            "smorf_index": smorf_index,
            "smorf_parent_map": smorf_parent_map,
            "hla_index": hla_index,
            "hla_evidence": hla_evidence,
        }

    def _args(self, inputs: dict[str, Path], resources: dict[str, Path], outdir: Path) -> argparse.Namespace:
        return argparse.Namespace(
            sample="CRYPTIC-T",
            cryptic_peptide_core=str(inputs["core"]),
            cryptic_peptide_parent_map=str(inputs["parent_map"]),
            upstream_manifest=str(inputs["manifest"]),
            outdir=str(outdir),
            resource_manifest=str(resources.get("manifest", "")),
            smorf_match_index=str(resources.get("smorf_index", "")),
            smorf_parent_map=str(resources.get("smorf_parent_map", "")),
            hla_ligand_match_index=str(resources.get("hla_index", "")),
            hla_ligand_evidence=str(resources.get("hla_evidence", "")),
            allow_missing_external_normal_resources=False,
            coordinate_matching_enabled=False,
        )

    def test_exact_sequence_and_class_matching(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            outdir = root / "08c"

            manifest = build_external_normal_qc(self._args(inputs, resources, outdir))

            self.assertEqual(manifest["external_normal_evaluation"], "evaluated")
            self.assertEqual(manifest["stage_counts"]["source_core_unique_peptides"], 5)
            self.assertEqual(manifest["stage_counts"]["annotated_core_unique_peptides"], 5)
            self.assertEqual(manifest["stage_counts"]["normal_smorf_exact_match_peptides"], 2)
            self.assertEqual(manifest["stage_counts"]["normal_hla_ligand_exact_match_peptides"], 2)
            self.assertEqual(manifest["stage_counts"]["both_external_normal_evidence_peptides"], 1)
            self.assertEqual(manifest["stage_counts"]["external_normal_excluded_unique_peptides"], 3)
            self.assertEqual(manifest["stage_counts"]["tumor_restricted_primary_core_unique_peptides"], 2)

            annotated = pd.read_csv(outdir / "cryptic_external_normal_annotated_core.tsv", sep="\t")
            retained = pd.read_csv(outdir / "cryptic_tumor_restricted_primary_core.tsv", sep="\t")
            excluded = pd.read_csv(outdir / "cryptic_external_normal_excluded.tsv", sep="\t")
            self.assertEqual(len(annotated), 5)
            self.assertEqual(len(retained) + len(excluded), len(annotated))
            self.assertEqual(set(retained["peptide_record_id"]), {"pep2", "pep5"})
            self.assertEqual(
                retained.set_index("peptide_record_id").loc["pep2", "normal_translation_status"],
                "normal_translation_not_detected_in_evaluated_smorf_atlas",
            )
            self.assertEqual(set(excluded["peptide_record_id"]), {"pep1", "pep3", "pep4"})
            both = excluded.set_index("peptide_record_id").loc["pep4", "external_normal_qc_reasons"]
            self.assertIn("normal_smorf_exact_peptide_match", both)
            self.assertIn("normal_hla_ligand_exact_peptide_match", both)

            fasta = (outdir / "cryptic_tumor_restricted_primary_core.fasta").read_text()
            self.assertIn(">pep2\nCDEFGHIK\n", fasta)
            self.assertIn(">pep5\nQRSTVWYA\n", fasta)
            self.assertNotIn("ACDEFGHI", fasta)

            evidence = pd.read_csv(outdir / "cryptic_external_normal_evidence.tsv", sep="\t")
            self.assertIn("normal_smorf_translation", set(evidence["evidence_type"]))
            self.assertIn("normal_hla_ligand", set(evidence["evidence_type"]))
            self.assertIn("source_core_parent_map", set(evidence["evidence_type"]))
            hla_evidence = evidence[evidence["evidence_type"].eq("normal_hla_ligand")]
            self.assertIn("HLA-A*02:01", set(hla_evidence["hla_ligand_donor_hla_genotype"]))
            self.assertTrue(
                evidence["normal_smorf_coordinate_status"]
                .eq("not_evaluable_coordinates_not_materialized")
                .all()
            )

            reused = build_external_normal_qc(self._args(inputs, resources, outdir))
            self.assertTrue(reused["reused"])

    def test_coordinate_matching_enabled_fails_in_v1(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = self._args(self._write_08b_core(root), self._write_resources(root), root / "08c")
            args.coordinate_matching_enabled = True
            with self.assertRaisesRegex(ValueError, "coordinate matching is not implemented in policy v1.0"):
                build_external_normal_qc(args)

    def test_upstream_manifest_missing_or_incomplete_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            inputs["manifest"].unlink()
            with self.assertRaisesRegex(FileNotFoundError, "upstream"):
                build_external_normal_qc(self._args(inputs, resources, root / "missing"))

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root, run_status="running")
            resources = self._write_resources(root)
            with self.assertRaisesRegex(ValueError, "run_status=complete"):
                build_external_normal_qc(self._args(inputs, resources, root / "incomplete"))

    def test_upstream_output_hash_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            inputs["core"].write_text(inputs["core"].read_text() + "CRYPTIC-T\tpepX\tMHC-I\t8\tACDEFGHI\n")

            with self.assertRaisesRegex(ValueError, "output_signature mismatch"):
                build_external_normal_qc(self._args(inputs, resources, root / "08c"))

    def test_parent_map_orphan_and_core_missing_parent_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core_rows = self.default_core_rows()
            parent_rows = self.default_parent_rows(core_rows)
            parent_rows.append(
                {
                    "sample": "CRYPTIC-T",
                    "peptide_record_id": "orphan",
                    "parent_record_id": "parent_orphan",
                    "source_parent_id": "source_orphan",
                    "mhc_class": "MHC-I",
                    "peptide_length": 8,
                    "peptide": "ACDEFGHI",
                }
            )
            inputs = self._write_08b_core(root, core_rows=core_rows, parent_rows=parent_rows)
            resources = self._write_resources(root)
            with self.assertRaisesRegex(ValueError, "orphan peptide_record_id"):
                build_external_normal_qc(self._args(inputs, resources, root / "orphan"))

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core_rows = self.default_core_rows()
            parent_rows = [row for row in self.default_parent_rows(core_rows) if row["peptide_record_id"] != "pep5"]
            inputs = self._write_08b_core(root, core_rows=core_rows, parent_rows=parent_rows)
            resources = self._write_resources(root)
            with self.assertRaisesRegex(ValueError, "missing parent-map"):
                build_external_normal_qc(self._args(inputs, resources, root / "missing_parent"))

    def test_noncore_parent_map_rows_with_blank_ids_are_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            parent_map = pd.read_csv(inputs["parent_map"], sep="\t", dtype=str, keep_default_na=False)
            parent_map["peptide_core_status"] = "core"
            parent_map.loc[len(parent_map)] = {
                "sample": "CRYPTIC-T",
                "peptide_record_id": "",
                "parent_record_id": "excluded_parent",
                "source_parent_id": "excluded_source",
                "mhc_class": "MHC-I",
                "peptide_length": "8",
                "peptide": "ACDEFGHI",
                "peptide_core_status": "excluded",
            }
            parent_map.to_csv(inputs["parent_map"], sep="\t", index=False)
            self.write_08b_manifest(inputs["manifest"], inputs["core"], inputs["parent_map"])

            manifest = build_external_normal_qc(self._args(inputs, resources, root / "08c"))

            evidence = pd.read_csv(root / "08c" / "cryptic_external_normal_evidence.tsv", sep="\t")
            self.assertEqual(manifest["stage_counts"]["source_core_unique_peptides"], 5)
            self.assertNotIn("excluded_parent", set(evidence["resource_record_id"]))

    def test_parent_map_sequence_class_length_and_sample_mismatches_fail(self) -> None:
        mutations = [
            ("sample", "OTHER", "sample mismatch"),
            ("peptide", "ACDEFGHK", "disagree"),
            ("peptide_length", 9, "peptide_length does not match"),
            ("mhc_class", "MHC-II", "HLA-II peptide length"),
        ]
        for column, value, pattern in mutations:
            with self.subTest(column=column), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                core_rows = self.default_core_rows()
                parent_rows = self.default_parent_rows(core_rows)
                parent_rows[0][column] = value
                inputs = self._write_08b_core(root, core_rows=core_rows, parent_rows=parent_rows)
                resources = self._write_resources(root)
                with self.assertRaisesRegex(ValueError, pattern):
                    build_external_normal_qc(self._args(inputs, resources, root / f"bad_{column}"))

    def test_invalid_sequence_class_length_and_sample_fail(self) -> None:
        mutations = [
            ("sample", "OTHER", "sample mismatch"),
            ("peptide", "ACDEFGHX", "non-standard"),
            ("peptide_length", 9, "peptide_length does not match"),
            ("mhc_class", "MHC-III", "unsupported mhc_class"),
        ]
        for column, value, pattern in mutations:
            with self.subTest(column=column), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                core_rows = self.default_core_rows()
                core_rows[0][column] = value
                parent_rows = self.default_parent_rows(core_rows)
                inputs = self._write_08b_core(root, core_rows=core_rows, parent_rows=parent_rows)
                resources = self._write_resources(root)
                with self.assertRaisesRegex(ValueError, pattern):
                    build_external_normal_qc(self._args(inputs, resources, root / f"bad_core_{column}"))

    def test_wrong_resource_contract_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root, freeze_version="wrong")
            with self.assertRaisesRegex(ValueError, "freeze_version"):
                build_external_normal_qc(self._args(inputs, resources, root / "wrong_freeze"))

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root, hla_i_lengths=[8, 9, 10])
            with self.assertRaisesRegex(ValueError, "HLA-I length policy"):
                build_external_normal_qc(self._args(inputs, resources, root / "wrong_hla_policy"))

    def test_resource_hash_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            with gzip.open(resources["smorf_index"], "at") as handle:
                handle.write("HLA-I\tMISMATCH\t8\n")

            with self.assertRaises(ValueError):
                build_external_normal_qc(self._args(inputs, resources, root / "08c"))

    def test_missing_resources_exploratory_output_is_not_primary_core(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = {"manifest": root / "missing_manifest.json"}
            args = self._args(inputs, resources, root / "08c_exploratory")
            args.allow_missing_external_normal_resources = True

            manifest = build_external_normal_qc(args)

            self.assertEqual(manifest["external_normal_evaluation"], "not_evaluated")
            self.assertEqual(manifest["binding_input_fasta"], "")
            self.assertEqual(manifest["binding_input_mode"], "not_available")
            self.assertEqual(manifest["stage_counts"]["tumor_restricted_primary_core_unique_peptides"], 0)
            annotated = pd.read_csv(args.outdir + "/cryptic_external_normal_annotated_core.tsv", sep="\t")
            primary = pd.read_csv(args.outdir + "/cryptic_tumor_restricted_primary_core.tsv", sep="\t")
            excluded = pd.read_csv(args.outdir + "/cryptic_external_normal_excluded.tsv", sep="\t")
            self.assertEqual(len(annotated), 5)
            self.assertEqual(len(primary), 0)
            self.assertEqual(len(excluded), len(annotated))
            self.assertTrue(excluded["external_normal_status"].eq("external_normal_resources_not_evaluated").all())

    def test_empty_core_can_complete_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root, core_rows=[], parent_rows=[])
            resources = self._write_resources(root)
            outdir = root / "empty"

            manifest = build_external_normal_qc(self._args(inputs, resources, outdir))

            self.assertEqual(manifest["stage_counts"]["source_core_unique_peptides"], 0)
            self.assertEqual(manifest["stage_counts"]["tumor_restricted_primary_core_unique_peptides"], 0)
            self.assertEqual((outdir / "cryptic_tumor_restricted_primary_core.fasta").read_text(), "")
            reused = build_external_normal_qc(self._args(inputs, resources, outdir))
            self.assertTrue(reused["reused"])

    def test_code_signature_change_invalidates_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._write_08b_core(root)
            resources = self._write_resources(root)
            outdir = root / "08c"
            args = self._args(inputs, resources, outdir)

            manifest = build_external_normal_qc(args)
            self.assertIn("script_sha256", manifest["input_signature"])

            manifest_path = outdir / "run_manifest.json"
            saved = json.loads(manifest_path.read_text())
            saved["input_signature"]["script_sha256"] = "different-script"
            manifest_path.write_text(json.dumps(saved, indent=2))

            with self.assertRaises(ValueError):
                build_external_normal_qc(args)


if __name__ == "__main__":
    unittest.main()
