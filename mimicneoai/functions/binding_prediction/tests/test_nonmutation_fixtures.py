from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from mimicneoai.functions.binding_prediction.nonmutation_workflow import (
    PVACBIND_COMPAT_COLUMNS,
    build_pvacbind_row,
    main,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "nonmutation"


class NonmutationAntigenFixtureTest(unittest.TestCase):
    def test_output_profiles_remain_fixed(self) -> None:
        self.assertEqual(len(PVACBIND_COMPAT_COLUMNS), 45)
        self.assertEqual(len(PVACBIND_COMPAT_COLUMNS), len(set(PVACBIND_COMPAT_COLUMNS)))

    def test_el_predictions_are_percentile_only_in_summary(self) -> None:
        epitope = {
            "Mutation": "source",
            "Sub-peptide Position": "1",
            "Epitope Seq": "ACDEFGHIK",
            "Peptide Length": "9",
            "mhc_class": "MHC-I",
        }
        predictions = {
            ("ACDEFGHIK", "HLA-A*02:01", "MHCflurry", "MHC-I", 9): {
                "ic50": "100",
                "percentile": "2",
            },
            ("ACDEFGHIK", "HLA-A*02:01", "MHCflurryEL", "MHC-I", 9): {
                "ic50": "1",
                "score": "0.9",
                "percentile": "0.1",
            },
        }
        row, _ = build_pvacbind_row(
            epitope,
            "HLA-A*02:01",
            ["MHCflurry", "MHCflurryEL"],
            predictions,
        )
        self.assertEqual(row["Best IC50 Score Method"], "MHCflurry")
        self.assertEqual(row["Best IC50 Score"], "100")
        self.assertEqual(row["Best Percentile Method"], "MHCflurryEL")
        self.assertEqual(row["Best Percentile"], "0.1")
        self.assertEqual(row["MHCflurryEL Presentation Score"], "0.9")

    def run_fixture(self, fasta_name: str) -> tuple[dict[str, object], list[dict[str, str]], list[dict[str, str]]]:
        with tempfile.TemporaryDirectory() as tempdir:
            outdir = Path(tempdir) / "out"
            self.assertEqual(
                main(
                    [
                        "-s",
                        "FIXTURE",
                        "--pep-fasta",
                        str(FIXTURE_DIR / fasta_name),
                        "--hla-file",
                        str(FIXTURE_DIR / "hla_final.result.txt"),
                        "-o",
                        str(outdir),
                        "--mhc-i-lengths",
                        "8,9",
                        "--mhc-ii-lengths",
                        "15",
                        "--algorithms",
                        "MHCflurry,MHCnuggetsII",
                        "--skip-prediction",
                    ]
                ),
                0,
            )
            summary = json.loads((outdir / "FIXTURE.mimicneoai_binding.summary.json").read_text())
            task_dir = outdir / "mimicneoai_epitope_tasks"
            window_manifest = json.loads((task_dir / "epitope_windows.manifest.json").read_text())
            self.assertEqual(window_manifest["input_signature"]["schema_version"], 2)
            with (task_dir / "epitope_windows.tsv").open(newline="") as handle:
                windows = list(csv.DictReader(handle, delimiter="\t"))
            with (task_dir / "binding_tasks.tsv").open(newline="") as handle:
                tasks = list(csv.DictReader(handle, delimiter="\t"))
            return summary, windows, tasks

    def test_cryptic_short_orf_and_duplicate_sources(self) -> None:
        summary, windows, tasks = self.run_fixture("cryptic.pep")
        self.assertEqual(summary["fasta_records"], 3)
        self.assertFalse(any(row["Mutation"] == "cryptic_short_orf" for row in windows))
        self.assertEqual(
            {row["Mutation"] for row in windows},
            {"cryptic_orf_1", "cryptic_orf_1_duplicate_source"},
        )
        unique_task_keys = {
            (row["peptide"], row["hla_allele"], row["algorithm"], row["mhc_class"])
            for row in tasks
        }
        self.assertEqual(len(unique_task_keys), len(tasks))
        self.assertLess(len(tasks), len(windows))

    def test_microbial_multiple_proteins_preserve_sources_and_deduplicate_tasks(self) -> None:
        summary, windows, tasks = self.run_fixture("microbial.fasta")
        self.assertEqual(summary["fasta_records"], 3)
        self.assertEqual(
            {row["Mutation"] for row in windows},
            {
                "taxon_1|protein_1",
                "taxon_2|protein_2_duplicate_sequence",
                "taxon_3|protein_3",
            },
        )
        duplicate_source_windows = [
            row
            for row in windows
            if row["Mutation"] in {"taxon_1|protein_1", "taxon_2|protein_2_duplicate_sequence"}
        ]
        self.assertTrue(duplicate_source_windows)
        self.assertEqual(len(tasks), summary["estimated_binding_task_rows"])
        self.assertLess(len(tasks), len(windows))


if __name__ == "__main__":
    unittest.main()
