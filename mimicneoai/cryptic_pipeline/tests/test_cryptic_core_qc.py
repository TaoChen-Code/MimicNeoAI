from __future__ import annotations

import argparse
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from mimicneoai.cryptic_pipeline.scripts.cryptic_core_qc import build_core, load_human_reference_matches


class CrypticCoreQCTest(unittest.TestCase):
    def _write_inputs(self, root: Path) -> dict[str, Path]:
        annot = pd.DataFrame(
            [
                {
                    "Name": "ENST_NC.1",
                    "nc_class": "noncoding",
                    "is_aberrant": True,
                    "TPM_tumor": 10.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 6.0,
                },
                {
                    "Name": "ENST_COD.1",
                    "nc_class": "coding_or_other",
                    "is_aberrant": True,
                    "TPM_tumor": 10.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 6.0,
                },
                {
                    "Name": "TRINITY_X|NOVEL|.|.",
                    "nc_class": "novel",
                    "is_aberrant": True,
                    "TPM_tumor": 8.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 5.0,
                },
                {
                    "Name": "ENST_BAD.1",
                    "nc_class": "noncoding",
                    "is_aberrant": True,
                    "TPM_tumor": 9.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 5.5,
                },
                {
                    "Name": "ENST_RELAXED.1",
                    "nc_class": "noncoding",
                    "is_aberrant": True,
                    "TPM_tumor": 1.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 5.0,
                },
                {
                    "Name": "ENST_DUP.1",
                    "nc_class": "noncoding",
                    "is_aberrant": True,
                    "TPM_tumor": 11.0,
                    "TPM_ctrl": 0.0,
                    "log2FC": 6.5,
                },
            ]
        )
        annot_path = root / "sample.aberrant_noncoding.annot.csv"
        annot.to_csv(annot_path, index=False)

        orf_ids = [
            "ENST_NC.1.p1",
            "ENST_COD.1.p1",
            "TRINITY_X|NOVEL|.|..p1",
            "ENST_BAD.1.p1",
            "ENST_RELAXED.1.p1",
            "ENST_DUP.1.p1",
        ]
        orf_final = pd.DataFrame({"TranscriptID": orf_ids})
        orf_final_path = root / "orf_final.csv"
        orf_final.to_csv(orf_final_path, index=False)

        fasta_path = root / "orf_filtered.fa"
        fasta_path.write_text(
            ">ENST_NC.1.p1\nACDEFGHIKLMNP\n"
            ">ENST_COD.1.p1\nACDEFGHIKLMNP\n"
            ">TRINITY_X|NOVEL|.|..p1\nKLMNPQRSTVWY*\n"
            ">ENST_BAD.1.p1\nACDEF*GHIKLMNP\n"
            ">ENST_RELAXED.1.p1\nTVWYACDEFGHIK\n"
            ">ENST_DUP.1.p1\nACDEFGHIKLMNP\n"
        )

        human_path = root / "human.fa"
        human_path.write_text(">HUMAN1\nTTTACDEFGHITTT\n")

        return {
            "annot": annot_path,
            "orf_final": orf_final_path,
            "fasta": fasta_path,
            "human": human_path,
        }

    def _args(self, paths: dict[str, Path], outdir: Path) -> argparse.Namespace:
        return argparse.Namespace(
            sample="TEST",
            ae_seps_fasta=str(paths["fasta"]),
            aeseps_annotation=str(paths["annot"]),
            orf_filtered_fasta=str(paths["fasta"]),
            orf_final=str(paths["orf_final"]),
            outdir=str(outdir),
            human_proteome_fasta=str(paths["human"]),
            allow_missing_human_reference=False,
            matched_control_sample="CTRL",
            reference_genome_fasta="",
            reference_gtf="",
            reference_lnc_gtf="",
            strandedness="reverse",
            min_tpm_tumor=5.0,
            max_tpm_ctrl=0.5,
            min_log2fc=4.0,
            mhc_i_lengths="8",
            mhc_ii_lengths="13",
        )

    def test_core_filters_and_resume_signature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            outdir = root / "core"

            manifest = build_core(self._args(paths, outdir))
            self.assertEqual(manifest["run_status"], "complete")

            parents = pd.read_csv(outdir / "cryptic_parent_records.tsv", sep="\t")
            core_parents = pd.read_csv(outdir / "cryptic_parent_core.tsv", sep="\t")
            excluded = pd.read_csv(outdir / "cryptic_parent_excluded.tsv", sep="\t")
            peptide_map = pd.read_csv(outdir / "cryptic_peptide_parent_map.tsv", sep="\t")

            self.assertEqual(
                set(core_parents["parent_record_id"]),
                {"ENST_NC.1.p1", "TRINITY_X|NOVEL|.|..p1", "ENST_DUP.1.p1"},
            )
            self.assertIn("ENST_COD.1.p1", set(excluded["parent_record_id"]))
            self.assertIn("ENST_BAD.1.p1", set(excluded["parent_record_id"]))
            self.assertIn("ENST_RELAXED.1.p1", set(excluded["parent_record_id"]))
            bad_reason = parents.set_index("parent_record_id").loc["ENST_BAD.1.p1", "parent_qc_reasons"]
            self.assertIn("internal_stop", bad_reason)
            relaxed_reason = parents.set_index("parent_record_id").loc["ENST_RELAXED.1.p1", "parent_qc_reasons"]
            self.assertIn("tumor_tpm_below_threshold", relaxed_reason)
            relaxed = parents.set_index("parent_record_id").loc["ENST_RELAXED.1.p1"]
            self.assertEqual(relaxed["expression_qc_status"], "fail")
            self.assertFalse(bool(relaxed["upstream_is_aberrant_consistent"]))

            matched = peptide_map[peptide_map["human_reference_peptide_status"].eq("human_reference_peptide_match")]
            self.assertEqual(set(matched["peptide"]), {"ACDEFGHI"})
            hla_i_fasta = (outdir / "cryptic_peptide_core_hla_i.fasta").read_text()
            self.assertNotIn("ACDEFGHI", hla_i_fasta)
            self.assertIn("KLMNPQRS", hla_i_fasta)
            hla_i_headers = [line for line in hla_i_fasta.splitlines() if line.startswith(">")]
            self.assertEqual(len(hla_i_headers), len(set(hla_i_headers)))
            self.assertTrue(all(header.startswith(">cryptic_core_MHC-I_") for header in hla_i_headers))

            reused = build_core(self._args(paths, outdir))
            self.assertTrue(reused["reused"])

            mutated_args = self._args(paths, outdir)
            mutated_args.min_tpm_tumor = 1.0
            with self.assertRaises(ValueError):
                build_core(mutated_args)

    def test_code_signature_change_invalidates_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            outdir = root / "code_signature"
            args = self._args(paths, outdir)

            manifest = build_core(args)
            self.assertIn("script_sha256", manifest["input_signature"])

            manifest_path = outdir / "run_manifest.json"
            saved = json.loads(manifest_path.read_text())
            saved["input_signature"]["script_sha256"] = "different-script"
            manifest_path.write_text(json.dumps(saved, indent=2))

            with self.assertRaises(ValueError):
                build_core(args)

    def test_zero_candidate_core_can_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            outdir = root / "zero"
            args = self._args(paths, outdir)
            args.min_tpm_tumor = 1000.0

            manifest = build_core(args)
            self.assertEqual(manifest["stage_counts"]["parent_core_records"], 0)
            self.assertEqual(manifest["stage_counts"]["unique_peptide_core_records"], 0)
            self.assertEqual((outdir / "cryptic_peptide_core.fasta").read_text(), "")
            reused = build_core(args)
            self.assertTrue(reused["reused"])

    def test_human_reference_keeps_standard_windows_around_u(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human = root / "human_with_u.fa"
            human.write_text(">HUMAN_U\nUUUACDEFGHIKUUUKLMNPQRS\n")

            matches, summary = load_human_reference_matches(
                human,
                {"ACDEFGHI", "KLMNPQRS", "FGHIKKLM"},
            )

            self.assertEqual(matches, {"ACDEFGHI", "KLMNPQRS"})
            self.assertEqual(summary["records_encountered"], 1)
            self.assertEqual(summary["records_with_noncanonical_residues"], 1)
            self.assertEqual(summary["standard_windows_evaluated"], 3)


if __name__ == "__main__":
    unittest.main()
