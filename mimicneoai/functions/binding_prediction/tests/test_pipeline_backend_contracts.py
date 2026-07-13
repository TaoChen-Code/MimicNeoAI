from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml

from mimicneoai.cryptic_pipeline import cryptic
from mimicneoai.microbial_pipeline.scripts import microbial_peptides
from mimicneoai.mutation_derived_pipeline import mutation_derived


REPO_ROOT = Path(__file__).resolve().parents[4]
CONFIG_DIR = REPO_ROOT / "mimicneoai" / "configures"


class PipelineBackendContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def test_packaged_configs_keep_pvactools_as_default(self) -> None:
        for filename in (
            "mutation_derived_configure.yaml",
            "cryptic_configure.yaml",
            "microbial_configure.yaml",
        ):
            with self.subTest(filename=filename):
                config = yaml.safe_load((CONFIG_DIR / filename).read_text())
                self.assertEqual(config["others"]["binding_prediction_backend"], "pvactools")
                self.assertIn("binding_prediction_algorithms", config["others"])
        mutation = yaml.safe_load((CONFIG_DIR / "mutation_derived_configure.yaml").read_text())
        self.assertNotIn("binding_prediction_max_task_rows", mutation["others"])
        self.assertEqual(
            mutation["others"]["binding_prediction_step_name"],
            "07.binding_prediction_mimicneoai",
        )
        for filename in ("cryptic_configure.yaml", "microbial_configure.yaml"):
            config = yaml.safe_load((CONFIG_DIR / filename).read_text())
            self.assertEqual(config["others"]["binding_prediction_max_task_rows"], 5_000_000)
            self.assertFalse(config["others"]["binding_prediction_force_large_samples"])

    def test_mutation_default_and_native_dispatch(self) -> None:
        config = {
            "path": {"output_dir": str(self.root)},
            "args": {"hla_binding_threads": 2},
            "step_name": {"annotation": "05.annotation", "hla": "06.hlatyping"},
            "others": {
                "QC": False,
                "hlatyping": False,
                "peptides_identification_and_binding_prediction": True,
                "tumor_with_matched_normal": True,
            },
        }
        paths = {"path": {"common": {"PVACTOOLS": "/tools/pvactools.sif"}}}
        tool = MagicMock()
        legacy_runner = MagicMock()
        with (
            patch.object(mutation_derived, "_variants_calling_and_annotation"),
            patch.object(mutation_derived, "Pvacseq", return_value=legacy_runner),
        ):
            mutation_derived._start_one_sample("TUMOR,NORMAL", config, paths, tool)
        legacy_runner.run_pvacseq_parallel.assert_called_once()
        tool.exec_cmd.assert_not_called()

        config["others"].update(
            {
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_step_name": "07.binding_prediction_mimicneoai_test",
                "binding_prediction_workers": 3,
            }
        )
        tool.reset_mock()
        with patch.object(mutation_derived, "_variants_calling_and_annotation"):
            mutation_derived._start_one_sample("TUMOR,NORMAL", config, paths, tool)
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("run_mimicneoai_binding_prediction.py", command)
        self.assertIn("07.binding_prediction_mimicneoai_test", command)
        self.assertIn("--workers 3", command)

        config["others"]["binding_prediction_step_name"] = "../invalid"
        with (
            patch.object(mutation_derived, "_variants_calling_and_annotation"),
            self.assertRaisesRegex(ValueError, "single directory name"),
        ):
            mutation_derived._start_one_sample("TUMOR,NORMAL", config, paths, tool)

    def test_mutation_runtime_tmp_directory_is_created(self) -> None:
        tmp_dir = self.root / "nested" / "mutation-tmp"
        config = {"path": {"tmp_dir": str(tmp_dir)}}

        mutation_derived._prepare_runtime_directories(config)

        self.assertTrue(tmp_dir.is_dir())
        self.assertEqual(config["path"]["tmp_dir"], str(tmp_dir.resolve()))

        with self.assertRaisesRegex(ValueError, "tmp_dir must be set"):
            mutation_derived._prepare_runtime_directories({"path": {"tmp_dir": ""}})

    def test_cryptic_default_and_native_dispatch(self) -> None:
        config = self.cryptic_config()
        paths = self.cryptic_paths()
        tool = MagicMock()
        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T", config, paths, tool)
        legacy_command = run_cmd.call_args.args[2]
        self.assertTrue(legacy_command[1].endswith("07-hla_binding_pred.py"))
        self.assertTrue(any(str(value).endswith("/07-hla_binding_pred") for value in legacy_command))

        config["others"].update(
            {
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_max_task_rows": 1234,
                "binding_prediction_force_large_samples": True,
            }
        )
        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T", config, paths, tool)
        native_command = run_cmd.call_args.args[2]
        self.assertTrue(native_command[1].endswith("07-hla_binding_pred_mimicneoai.py"))
        self.assertIn("--max-task-rows", native_command)
        self.assertIn("1234", native_command)
        self.assertIn("--force-large-samples", native_command)

    def test_microbial_default_and_native_dispatch(self) -> None:
        sample = "MICROBIAL-T"
        peptide_dir = self.root / sample / "06.MicrobialPeptidesIdentification"
        peptide_dir.mkdir(parents=True)
        (peptide_dir / f"{sample}.peptide.fasta").write_text(">protein\nACDEFGHIKLMNPQR\n")
        config = {
            "path": {"output_dir": str(self.root)},
            "args": {"hla_binding_threads": 2},
            "step_name": {
                "blastx": "06.MicrobialPeptidesIdentification",
                "pvacbind": "08.MicrobialPeptidesBindingPrediction",
                "hla": "07.HLA-HD",
            },
            "others": {},
        }
        tool = MagicMock()
        with patch.object(microbial_peptides, "pvacbind") as legacy:
            microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, {}, tool)
        legacy.assert_called_once()
        tool.exec_cmd.assert_not_called()

        config["others"].update(
            {
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_max_task_rows": 4321,
                "binding_prediction_force_large_samples": True,
            }
        )
        tool.reset_mock()
        microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, {}, tool)
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("hla_binding_pred_mimicneoai.py", command)
        self.assertIn("08.MicrobialPeptidesBindingPrediction_mimicneoai", command)
        self.assertIn("--max-task-rows 4321", command)
        self.assertIn("--force-large-samples", command)

    def cryptic_config(self) -> dict[str, object]:
        return {
            "path": {"input_dir": str(self.root / "input"), "output_dir": str(self.root)},
            "args": {"threads": 2, "hla_binding_threads": 2},
            "others": {
                "QC": False,
                "alignment": False,
                "known": False,
                "novel": False,
                "salmon_quant": False,
                "salmon_quant_control": False,
                "hlatyping": False,
                "extract_aeseps": False,
                "hla_binding_pred": True,
            },
        }

    @staticmethod
    def cryptic_paths() -> dict[str, object]:
        return {
            "path": {
                "cryptic": {"TRINITY_SIF": "/tools/trinity.sif"},
                "common": {"PVACTOOLS": "/tools/pvactools.sif"},
            },
            "database": {
                "cryptic": {
                    "STAR_GENOME_DIR": "/ref/star",
                    "REF": {
                        "REF_DIR": "/ref",
                        "REF_FA": "/ref/genome.fa",
                        "REF_GTF": "/ref/gencode.gtf",
                        "REF_LNC_GTF": "/ref/lnc.gtf",
                    },
                },
                "common": {
                    "HLA": {
                        "FREQ_DATA_DIR": "/ref/hla/freq",
                        "HLA_GENE": "/ref/hla/gene",
                        "DICTIONARY": "/ref/hla/dictionary",
                        "BOWTIE2_INDEX": "/ref/hla/index",
                    }
                },
            },
        }


if __name__ == "__main__":
    unittest.main()
