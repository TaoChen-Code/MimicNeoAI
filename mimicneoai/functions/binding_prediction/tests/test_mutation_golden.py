from __future__ import annotations

import csv
import importlib
import json
import tempfile
import unittest
from pathlib import Path
from typing import List, Optional

from mimicneoai.functions.binding_prediction.adapters.base import AdapterConfig, reusable_normalized_output
from mimicneoai.functions.binding_prediction.adapters.mhcnuggets import (
    MhcnuggetsAdapter,
    build_mhcnuggets_command,
    human_proteome_percentile,
    reference_rank_fraction,
)
from mimicneoai.functions.binding_prediction.allele_support import (
    AlleleSupportMatrix,
    SupportCatalog,
    allele_key,
    catalog_from_lines,
)
from mimicneoai.functions.binding_prediction.runner import merge_predictions, unsupported_reason
from mimicneoai.functions.binding_prediction.schema import (
    PREDICTION_FIELDS,
    BindingPrediction,
    BindingTask,
    PredictionJob,
    write_prediction_rows,
)
from mimicneoai.mutation_derived_pipeline.scripts.mutation_epitope_prediction.hla_parser import (
    parse_hlahd_result,
)
from mimicneoai.mutation_derived_pipeline.scripts.mutation_epitope_prediction.sequence_utils import (
    ProteinPair,
    generate_epitope_windows,
    locate_mutation_region,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "mutation_golden"
PREP_MODULE = importlib.import_module(
    "mimicneoai.mutation_derived_pipeline.scripts.mutation_epitope_prediction.00_prepare_pvacseq_sources"
)
TASK_MODULE = importlib.import_module(
    "mimicneoai.mutation_derived_pipeline.scripts.mutation_epitope_prediction.01_build_epitope_tasks"
)
MERGE_MODULE = importlib.import_module(
    "mimicneoai.mutation_derived_pipeline.scripts.mutation_epitope_prediction.02_merge_binding_predictions"
)


def write_tsv(
    path: Path,
    rows: List[dict[str, object]],
    fieldnames: Optional[List[str]] = None,
) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class MutationWindowGoldenTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixture = json.loads((FIXTURE_DIR / "events.json").read_text())

    def test_mutation_windows_match_golden_fixture(self) -> None:
        for event in self.fixture["events"]:
            with self.subTest(event=event["name"]):
                pair = ProteinPair(
                    event_id=event["name"],
                    wt_sequence=event["wt_sequence"],
                    mt_sequence=event["mt_sequence"],
                )
                self.assertEqual(
                    locate_mutation_region(pair.wt_sequence, pair.mt_sequence),
                    tuple(event["mutation_region"]),
                )
                windows = generate_epitope_windows(
                    pair,
                    epitope_lengths=(8, 9, 15),
                    extended_length=27,
                    variant_type=event["variant_type"],
                )
                for length_text, expected in event["window_summary"].items():
                    selected = [window for window in windows if window.length == int(length_text)]
                    observed = [len(selected), min(w.start for w in selected), max(w.start for w in selected)]
                    self.assertEqual(observed, expected)
                    self.assertTrue(all(w.peptide != w.wt_peptide for w in selected))

    def test_frameshift_wt_is_not_scheduled_for_prediction(self) -> None:
        rows = [
            {
                "event_id": "fs_event",
                "variant_type": "FS",
                "mhc_class": "MHC-I",
                "peptide_length": 9,
                "mt_epitope_seq": "LMNPQRSTV",
                "wt_epitope_seq": "LMNPQASTV",
            },
            {
                "event_id": "snv_event",
                "variant_type": "missense",
                "mhc_class": "MHC-I",
                "peptide_length": 9,
                "mt_epitope_seq": "ACDEFGHIK",
                "wt_epitope_seq": "ACDEYGHIK",
            },
        ]
        prediction_rows = TASK_MODULE.build_prediction_peptide_rows(rows)
        by_peptide = {row["peptide"]: row["source_types"] for row in prediction_rows}
        self.assertEqual(by_peptide["LMNPQRSTV"], "MT")
        self.assertNotIn("LMNPQASTV", by_peptide)
        self.assertEqual(by_peptide["ACDEFGHIK"], "MT")
        self.assertEqual(by_peptide["ACDEYGHIK"], "WT")

    def test_pvactools_source_manifest_rejects_changed_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            manifest = Path(tempdir) / "source.manifest.json"
            signature = {"sample": "TEST", "input_vcf": {"size": 10}}
            PREP_MODULE.write_source_manifest(manifest, signature)
            PREP_MODULE.validate_source_manifest(manifest, signature)
            with self.assertRaisesRegex(RuntimeError, "source inputs differ"):
                PREP_MODULE.validate_source_manifest(
                    manifest,
                    {"sample": "TEST", "input_vcf": {"size": 11}},
                )

    def test_step01_builds_golden_fixture_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            converter_path = root / "converter.tsv"
            fasta_path = root / "proteins.fasta"
            outdir = root / "tasks"
            variant_details = {
                "snv": ("missense", "N/A", "12"),
                "inframe_insertion": ("inframe_ins", "-/AAA", "12"),
                "inframe_deletion": ("inframe_del", "MNP/-", "11-13"),
                "frameshift": ("FS", "L/Y", "10"),
                "n_terminal_snv": ("missense", "A/V", "1"),
            }
            converter_rows = []
            fasta_lines = []
            for offset, event in enumerate(self.fixture["events"], start=1):
                index = f"EVENT.{offset}"
                variant_type, amino_acid_change, protein_position = variant_details[event["name"]]
                converter_rows.append(
                    {
                        "index": index,
                        "chromosome_name": "1",
                        "start": str(1000 + offset),
                        "stop": str(1000 + offset),
                        "reference": "A",
                        "variant": "T",
                        "gene_name": f"GENE{offset}",
                        "transcript_name": f"ENST_TEST_{offset}",
                        "hgvsc": f"c.{offset}A>T",
                        "hgvsp": f"p.TEST{offset}",
                        "variant_type": variant_type,
                        "protein_position": protein_position,
                        "amino_acid_change": amino_acid_change,
                    }
                )
                fasta_lines.extend(
                    [
                        f">WT.{index}",
                        event["wt_sequence"],
                        f">MT.{index}",
                        event["mt_sequence"],
                    ]
                )
            write_tsv(converter_path, converter_rows)
            fasta_path.write_text("\n".join(fasta_lines) + "\n")

            self.assertEqual(
                TASK_MODULE.main(
                    [
                        "-s",
                        "GOLDEN",
                        "--converter-tsv",
                        str(converter_path),
                        "--protein-fasta",
                        str(fasta_path),
                        "--hla-file",
                        str(FIXTURE_DIR / "hla_final.result.txt"),
                        "--mhc-i-lengths",
                        "8",
                        "--mhc-ii-lengths",
                        "15",
                        "--algorithms",
                        "MHCflurry,MHCnuggetsII",
                        "-o",
                        str(outdir),
                    ]
                ),
                0,
            )
            summary = json.loads((outdir / "summary.json").read_text())
            self.assertEqual(summary["variant_annotation_rows"], 5)
            self.assertEqual(summary["fasta_pair_status"], {"fasta_pair_found": 5})

            with (outdir / "epitope_windows.tsv").open(newline="") as handle:
                windows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["variant_type"] for row in windows},
                {"missense", "inframe_ins", "inframe_del", "FS"},
            )
            n_terminal = [row for row in windows if row["protein_position"] == "1"]
            self.assertTrue(n_terminal)
            self.assertTrue(all(row["window_start_0based"] == "0" for row in n_terminal))

            with (outdir / "prediction_peptides.tsv").open(newline="") as handle:
                prediction_peptides = list(csv.DictReader(handle, delimiter="\t"))
            fs_mt_peptides = {
                row["mt_epitope_seq"]
                for row in windows
                if row["variant_type"] == "FS"
            }
            fs_wt_peptides = {
                row["wt_epitope_seq"]
                for row in windows
                if row["variant_type"] == "FS" and row["wt_epitope_seq"]
            }
            non_fs_peptides = {
                peptide
                for row in windows
                if row["variant_type"] != "FS"
                for peptide in (row["mt_epitope_seq"], row["wt_epitope_seq"])
                if peptide
            }
            fs_only_wt_peptides = fs_wt_peptides - fs_mt_peptides - non_fs_peptides
            scheduled = {row["peptide"]: row["source_types"] for row in prediction_peptides}
            self.assertTrue(fs_mt_peptides.issubset(scheduled))
            self.assertTrue(fs_only_wt_peptides)
            self.assertTrue(fs_only_wt_peptides.isdisjoint(scheduled))

            with (outdir / "binding_tasks.tsv").open(newline="") as handle:
                tasks = list(csv.DictReader(handle, delimiter="\t"))
            class_ii_alleles = {row["hla_allele"] for row in tasks if row["mhc_class"] == "MHC-II"}
            self.assertIn("DPA1*01:03-DPB1*02:01", class_ii_alleles)
            self.assertIn("DQA1*05:01-DQB1*04:02", class_ii_alleles)


class HlaAndAdapterContractTest(unittest.TestCase):
    def test_hlahd_dq_dp_pairing_and_dra_exclusion(self) -> None:
        alleles = parse_hlahd_result(FIXTURE_DIR / "hla_final.result.txt")
        self.assertEqual(alleles.mhc_i, ("HLA-A*02:01", "HLA-A*11:01", "HLA-B*40:01"))
        self.assertEqual(
            set(alleles.mhc_ii),
            {
                "DRB1*04:05",
                "DRB1*11:01",
                "DPA1*01:03-DPB1*02:01",
                "DPA1*01:03-DPB1*04:01",
                "DQA1*03:01-DQB1*03:01",
                "DQA1*03:01-DQB1*04:02",
                "DQA1*05:01-DQB1*03:01",
                "DQA1*05:01-DQB1*04:02",
            },
        )
        self.assertFalse(any("DRA" in allele for allele in alleles.mhc_ii))

    def test_mhcnuggets_percentile_scale_and_command(self) -> None:
        self.assertEqual(human_proteome_percentile("0.8719"), "87.19")
        self.assertEqual(human_proteome_percentile(""), "")
        self.assertEqual(
            reference_rank_fraction(
                20.0,
                first_percentile=10.0,
                downsampled=[10.0, 20.0, 20.0, 30.0],
                first_values=[1.0, 5.0, 10.0],
                hp_length=100,
            ),
            0.625,
        )
        job = PredictionJob(
            algorithm="MHCnuggetsI",
            mhc_class="MHC-I",
            hla_allele="HLA-A*02:01",
            peptide_length=9,
            chunk_index=1,
            peptides=("ACDEFGHIK",),
            outdir="/tmp/not-used",
        )
        command = build_mhcnuggets_command(job, rank_output=True, config=AdapterConfig())
        self.assertNotIn("-M", command)
        self.assertEqual(command[-2:], ["-r", "True"])

    def test_unsupported_alleles_are_explicitly_skipped(self) -> None:
        support = AlleleSupportMatrix(
            AdapterConfig(),
            catalogs={
                "MHCnuggetsI": SupportCatalog(frozenset({allele_key("HLA-A*02:01")}), "fixture"),
                "NetMHCpan": SupportCatalog(frozenset({allele_key("HLA-A*02:01")}), "fixture"),
                "NNalign": SupportCatalog(frozenset({allele_key("DRB1*11:01")}), "fixture"),
            },
        )
        cases = [
            (BindingTask("ACDEFGHIK", "HLA-F*01:01", "MHCnuggetsI", "MHC-I"), True),
            (BindingTask("ACDEFGHIK", "HLA-A*11:353", "NetMHCpan", "MHC-I"), True),
            (BindingTask("ACDEFGHIKLMNPQR", "DQA1*05:05-DQB1*04:01", "NNalign", "MHC-II"), True),
            (BindingTask("ACDEFGHIK", "HLA-A*02:01", "MHCnuggetsI", "MHC-I"), False),
            (BindingTask("ACDEFGHIKLMNPQR", "DRB1*11:01", "NNalign", "MHC-II"), False),
        ]
        for task, expected_skipped in cases:
            with self.subTest(task=task):
                self.assertEqual(bool(unsupported_reason(task, support)), expected_skipped)

    def test_skipped_task_is_retained_in_normalized_output(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "predictions.tsv"
            task = BindingTask("ACDEFGHIK", "HLA-A*99:99", "NetMHCpan", "MHC-I")
            summary = merge_predictions(
                [],
                output,
                [(task, "unsupported_allele_by_predictor")],
            )
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["status_counts"], {"skipped": 1})
            self.assertEqual(rows[0]["status"], "skipped")
            self.assertEqual(rows[0]["error"], "unsupported_allele_by_predictor")

    def test_predictor_allele_spellings_share_one_support_key(self) -> None:
        self.assertEqual(allele_key("HLA-A*02:01"), allele_key("HLA-A0201"))
        self.assertEqual(allele_key("DRB1*11:01"), allele_key("DRB1_1101"))
        self.assertEqual(
            allele_key("DPA1*01:03-DPB1*02:01"),
            allele_key("HLA-DPA10103-DPB10201"),
        )

    def test_support_catalog_parsing_and_unknown_catalog_fail_open(self) -> None:
        catalog = catalog_from_lines(
            ["# supported", "HLA-A02:01", "HLA-B*44:03", "", "----"],
            "fixture",
        )
        self.assertEqual(
            catalog.alleles,
            frozenset({allele_key("HLA-A*02:01"), allele_key("HLA-B*44:03")}),
        )
        support = AlleleSupportMatrix(
            AdapterConfig(),
            catalogs={"NetMHCpan": SupportCatalog(None, "fixture", error="probe failed")},
        )
        task = BindingTask("ACDEFGHIK", "HLA-A*99:99", "NetMHCpan", "MHC-I")
        self.assertEqual(unsupported_reason(task, support), "")

    def test_mhcnuggets_class_i_length_limit_precedes_exact_model(self) -> None:
        support = AlleleSupportMatrix(
            AdapterConfig(),
            catalogs={
                "MHCnuggetsI": SupportCatalog(
                    frozenset({allele_key("HLA-A*02:01")}),
                    "fixture",
                )
            },
        )
        task = BindingTask("ACDEFGHIKLMNPQRST", "HLA-A*02:01", "MHCnuggetsI", "MHC-I")
        self.assertEqual(unsupported_reason(task, support), "unsupported_allele_by_predictor")

    def test_mhcnuggets_resume_reuses_normalized_chunk(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            job = PredictionJob(
                algorithm="MHCnuggetsI",
                mhc_class="MHC-I",
                hla_allele="HLA-A*02:01",
                peptide_length=9,
                chunk_index=1,
                peptides=("ACDEFGHIK",),
                outdir=tempdir,
            )
            job.normalized_path.parent.mkdir(parents=True)
            write_prediction_rows(
                job.normalized_path,
                [
                    BindingPrediction(
                        peptide="ACDEFGHIK",
                        hla_allele="HLA-A*02:01",
                        algorithm="MHCnuggetsI",
                        mhc_class="MHC-I",
                        peptide_length=9,
                        ic50="100",
                    )
                ],
            )
            sentinel = job.normalized_path.read_text()
            adapter = MhcnuggetsAdapter(
                AdapterConfig(
                    resume=True,
                    mhcnuggets_python_bin="/path/that/must/not/run",
                    mhcnuggets_script="/path/that/must/not/run",
                )
            )
            self.assertEqual(adapter.run_job(job), job.normalized_path)
            self.assertEqual(job.normalized_path.read_text(), sentinel)

            changed_job = PredictionJob(
                algorithm="MHCnuggetsI",
                mhc_class="MHC-I",
                hla_allele="HLA-A*02:01",
                peptide_length=9,
                chunk_index=1,
                peptides=("LMNPQRSTV",),
                outdir=tempdir,
            )
            self.assertFalse(reusable_normalized_output(changed_job, resume=True))


class MergeGoldenTest(unittest.TestCase):
    def test_pvacseq_output_profile_remains_fixed(self) -> None:
        self.assertEqual(len(MERGE_MODULE.PVACSEQ_COLUMNS), 111)
        self.assertEqual(len(MERGE_MODULE.PVACSEQ_COLUMNS), len(set(MERGE_MODULE.PVACSEQ_COLUMNS)))

    def test_el_predictions_do_not_enter_ic50_or_fold_summaries(self) -> None:
        epitope = {
            "event_id": "event",
            "pvacseq_index": "SNV.EL",
            "variant_type": "missense",
            "mhc_class": "MHC-I",
            "peptide_length": "9",
            "window_start_0based": "0",
            "window_end_0based": "9",
            "mutation_start_0based": "4",
            "mutation_end_0based": "5",
            "mt_epitope_seq": "ACDEFGHIK",
            "wt_epitope_seq": "ACDEYGHIK",
        }
        annotation = {"amino_acid_change": "Y/F", "variant_type": "missense"}
        predictions = {
            ("ACDEFGHIK", "HLA-A*02:01", "MHCflurry", "MHC-I", 9): {
                "ic50": "100",
                "percentile": "2",
                "status": "ok",
            },
            ("ACDEYGHIK", "HLA-A*02:01", "MHCflurry", "MHC-I", 9): {
                "ic50": "400",
                "percentile": "4",
                "status": "ok",
            },
            ("ACDEFGHIK", "HLA-A*02:01", "MHCflurryEL", "MHC-I", 9): {
                "ic50": "1",
                "score": "0.9",
                "percentile": "0.1",
                "status": "ok",
            },
            ("ACDEYGHIK", "HLA-A*02:01", "MHCflurryEL", "MHC-I", 9): {
                "ic50": "2",
                "score": "0.8",
                "percentile": "0.2",
                "status": "ok",
            },
        }
        algorithms = ["MHCflurry", "MHCflurryEL"]
        output_fields = MERGE_MODULE.PVACSEQ_COLUMNS
        row, _ = MERGE_MODULE.build_output_row(
            epitope,
            annotation,
            "HLA-A*02:01",
            algorithms,
            predictions,
            {("SNV.EL", "MHC-I", 9): 0},
            output_fields,
        )
        self.assertEqual(row["Best MT IC50 Score Method"], "MHCflurry")
        self.assertEqual(row["Best MT IC50 Score"], "100")
        self.assertEqual(row["Corresponding Fold Change"], "4")
        self.assertEqual(row["Best MT Percentile Method"], "MHCflurryEL")
        self.assertEqual(row["Best MT Percentile"], "0.1")
        self.assertEqual(row["MHCflurryEL Presentation MT Score"], "0.9")
        self.assertNotIn("MHCflurryEL Fold Change", row)

    def test_best_median_fold_change_and_fs_wt_blank(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            annotation_path = root / "variant_annotation.tsv"
            windows_path = root / "epitope_windows.tsv"
            tasks_path = root / "binding_tasks.tsv"
            predictions_path = root / "binding_predictions.long.tsv"
            output_path = root / "merged.tsv"

            write_tsv(
                annotation_path,
                [
                    {"pvacseq_index": "SNV.1", "variant_type": "missense", "gene_name": "GENE1"},
                    {"pvacseq_index": "FS.1", "variant_type": "FS", "gene_name": "GENE2"},
                ],
            )
            common_window = {
                "mhc_class": "MHC-I",
                "peptide_length": 9,
                "window_start_0based": 0,
                "window_end_0based": 9,
                "mutation_start_0based": 4,
                "mutation_end_0based": 5,
                "extended_mt_epitope_seq": "ACDEFGHIKLMNPQRSTVWYACDEF",
                "extended_length": 25,
                "transcript_name": "ENST_TEST",
                "amino_acid_change": "Y/F",
                "protein_position": "5",
            }
            write_tsv(
                windows_path,
                [
                    {
                        **common_window,
                        "event_id": "snv_event",
                        "pvacseq_index": "SNV.1",
                        "variant_type": "missense",
                        "mt_epitope_seq": "ACDEFGHIK",
                        "wt_epitope_seq": "ACDEYGHIK",
                    },
                    {
                        **common_window,
                        "event_id": "fs_event",
                        "pvacseq_index": "FS.1",
                        "variant_type": "FS",
                        "mt_epitope_seq": "LMNPQRSTV",
                        "wt_epitope_seq": "LMNPQASTV",
                    },
                ],
            )
            task_rows = []
            for peptide in ("ACDEFGHIK", "ACDEYGHIK", "LMNPQRSTV"):
                for algorithm in ("MHCflurry", "NetMHCpan"):
                    task_rows.append(
                        {
                            "peptide": peptide,
                            "hla_allele": "HLA-A*02:01",
                            "algorithm": algorithm,
                            "mhc_class": "MHC-I",
                        }
                    )
            write_tsv(tasks_path, task_rows)

            values = {
                ("ACDEFGHIK", "MHCflurry"): (100, 1),
                ("ACDEYGHIK", "MHCflurry"): (400, 4),
                ("ACDEFGHIK", "NetMHCpan"): (200, 2),
                ("ACDEYGHIK", "NetMHCpan"): (300, 3),
                ("LMNPQRSTV", "MHCflurry"): (80, 0.8),
                ("LMNPQRSTV", "NetMHCpan"): (120, 1.2),
            }
            prediction_rows = []
            for (peptide, algorithm), (ic50, percentile) in values.items():
                prediction_rows.append(
                    {
                        "peptide": peptide,
                        "hla_allele": "HLA-A*02:01",
                        "algorithm": algorithm,
                        "mhc_class": "MHC-I",
                        "peptide_length": 9,
                        "ic50": ic50,
                        "percentile": percentile,
                        "status": "ok",
                    }
                )
            write_tsv(predictions_path, prediction_rows, PREDICTION_FIELDS)

            self.assertEqual(
                MERGE_MODULE.main(
                    [
                        "--variant-annotation",
                        str(annotation_path),
                        "--epitope-windows",
                        str(windows_path),
                        "--binding-predictions",
                        str(predictions_path),
                        "--binding-tasks",
                        str(tasks_path),
                        "--output-profile",
                        "extended",
                        "-o",
                        str(output_path),
                    ]
                ),
                0,
            )
            with output_path.open(newline="") as handle:
                rows = {row["Index"]: row for row in csv.DictReader(handle, delimiter="\t")}

            snv = rows["SNV.1"]
            self.assertEqual(snv["Best MT IC50 Score Method"], "MHCflurry")
            self.assertEqual(snv["Best MT IC50 Score"], "100")
            self.assertEqual(snv["Corresponding WT IC50 Score"], "400")
            self.assertEqual(snv["Corresponding Fold Change"], "4")
            self.assertEqual(snv["Median MT IC50 Score"], "150")
            self.assertEqual(snv["Median WT IC50 Score"], "350")
            self.assertEqual(snv["Median Fold Change"], "2.75")
            self.assertEqual(snv["Best MT Percentile"], "1")
            self.assertEqual(snv["Corresponding WT Percentile"], "4")
            self.assertEqual(snv["Median MT Percentile"], "1.5")
            self.assertEqual(snv["Median WT Percentile"], "3.5")

            frameshift = rows["FS.1"]
            self.assertEqual(frameshift["WT Epitope Seq"], "")
            self.assertEqual(frameshift["Corresponding WT IC50 Score"], "")
            self.assertEqual(frameshift["Corresponding Fold Change"], "")
            self.assertEqual(frameshift["Median WT IC50 Score"], "")
            self.assertEqual(frameshift["Best MT IC50 Score"], "80")
            self.assertEqual(frameshift["Median MT IC50 Score"], "100")


if __name__ == "__main__":
    unittest.main()
