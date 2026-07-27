from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml

from mimicneoai.cryptic_pipeline import cryptic
from mimicneoai.microbial_pipeline import microbial
from mimicneoai.microbial_pipeline.scripts import microbial_peptides
from mimicneoai.mutation_derived_pipeline import mutation_derived


REPO_ROOT = Path(__file__).resolve().parents[4]
CONFIG_DIR = REPO_ROOT / "mimicneoai" / "configures"


class PipelineBackendContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def test_packaged_configs_set_expected_binding_defaults(self) -> None:
        mutation = yaml.safe_load((CONFIG_DIR / "mutation_derived_configure.yaml").read_text())
        self.assertEqual(mutation["others"]["binding_prediction_backend"], "mimicneoai")
        self.assertEqual(mutation["others"]["binding_prediction_preset"], "fast")
        self.assertEqual(mutation["others"]["mhc_i_epitope_lengths"], "8,9,10,11")
        self.assertEqual(mutation["others"]["mhc_ii_epitope_lengths"], "13,14,15,16,17")
        self.assertNotIn("NNalign", mutation["others"]["binding_prediction_algorithms"])
        self.assertNotIn("binding_prediction_max_task_rows", mutation["others"])
        self.assertEqual(
            mutation["others"]["binding_prediction_step_name"],
            "07.binding_prediction_mimicneoai",
        )

        cryptic = yaml.safe_load((CONFIG_DIR / "cryptic_configure.yaml").read_text())
        self.assertTrue(cryptic["others"]["cryptic_core_qc"])
        self.assertEqual(cryptic["others"]["binding_prediction_backend"], "mimicneoai")
        self.assertEqual(cryptic["others"]["binding_prediction_preset"], "fast")
        self.assertEqual(cryptic["others"]["mhcI_lengths"], "8,9,10,11")
        self.assertEqual(cryptic["others"]["mhcII_lengths"], "13,14,15,16,17")
        self.assertFalse(cryptic["others"]["allow_missing_human_reference"])
        self.assertIn("binding_prediction_algorithms", cryptic["others"])
        self.assertNotIn("NNalign", cryptic["others"]["binding_prediction_algorithms"])
        self.assertEqual(cryptic["others"]["binding_prediction_max_task_rows"], 5_000_000)
        self.assertFalse(cryptic["others"]["binding_prediction_force_large_samples"])

        microbial = yaml.safe_load((CONFIG_DIR / "microbial_configure.yaml").read_text())
        self.assertEqual(microbial["others"]["binding_prediction_backend"], "mimicneoai")
        self.assertEqual(microbial["others"]["binding_prediction_preset"], "fast")
        self.assertNotIn("NNalign", microbial["others"]["binding_prediction_algorithms"])
        self.assertEqual(
            microbial["others"]["binding_prediction_step_name"],
            "08.MicrobialPeptidesBindingPrediction_mimicneoai",
        )
        self.assertEqual(microbial["others"]["binding_prediction_max_task_rows"], 5_000_000)
        self.assertFalse(microbial["others"]["binding_prediction_force_large_samples"])
        self.assertFalse(microbial["others"]["tumor_with_matched_normal"])
        self.assertEqual(
            microbial["others"]["protein_hit_qc_policy_version"],
            "microbial_protein_hit_qc_v1.0",
        )

        paths = yaml.safe_load((CONFIG_DIR / "paths.yaml").read_text())
        common_paths = paths["path"]["common"]
        self.assertTrue(common_paths["APPTAINER_BIN"].endswith("/apptainer"))
        self.assertTrue(common_paths["BCFTOOLS_BIN"].endswith("/bcftools"))
        self.assertTrue(common_paths["TABIX_BIN"].endswith("/tabix"))
        predictor_paths = paths["path"]["common"]["BINDING_PREDICTORS"]
        self.assertTrue(
            predictor_paths["MHCFLURRY_PREDICT_BIN"].endswith("mhcflurry-predict")
        )
        self.assertTrue(
            predictor_paths["NETMHCPAN_BIN"].endswith("netMHCpan")
        )
        self.assertTrue(
            predictor_paths["IEDB_MHCII_SCRIPT"].endswith("mhc_II_binding.py")
        )

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
        legacy_runner.run_pvacseq_parallel.assert_not_called()
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("run_mimicneoai_binding_prediction.py", command)
        self.assertIn("07.binding_prediction_mimicneoai", command)
        self.assertIn("--preset fast", command)
        self.assertIn("--apptainer", command)
        self.assertIn("--bcftools", command)
        self.assertIn("--tabix", command)
        self.assertNotIn("NNalign", command)

        config["others"]["binding_prediction_backend"] = "pvactools"
        tool.reset_mock()
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
                "binding_prediction_preset": "fast",
                "binding_prediction_workers": 3,
            }
        )
        tool.reset_mock()
        with patch.object(mutation_derived, "_variants_calling_and_annotation"):
            mutation_derived._start_one_sample("TUMOR,NORMAL", config, paths, tool)
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("run_mimicneoai_binding_prediction.py", command)
        self.assertIn("07.binding_prediction_mimicneoai_test", command)
        self.assertIn("--preset fast", command)
        self.assertIn("--workers 3", command)

        paths["path"]["common"].update(
            {
                "APPTAINER_BIN": "/tools/apptainer",
                "BCFTOOLS_BIN": "/tools/bcftools",
                "TABIX_BIN": "/tools/tabix",
            }
        )
        paths["path"]["common"]["BINDING_PREDICTORS"] = {
            "MHCFLURRY_DOWNLOADS_DIR": "/tools/mhcflurry-models",
            "MHCNUGGETS_CWD": "/tools/mhcnuggets",
            "NETMHCPAN_BIN": "/tools/netMHCpan",
        }
        tool.reset_mock()
        with patch.object(mutation_derived, "_variants_calling_and_annotation"):
            mutation_derived._start_one_sample("TUMOR,NORMAL", config, paths, tool)
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn(
            "--mhcflurry-downloads-dir /tools/mhcflurry-models", command
        )
        self.assertIn("--apptainer /tools/apptainer", command)
        self.assertIn("--bcftools /tools/bcftools", command)
        self.assertIn("--tabix /tools/tabix", command)
        self.assertIn("--mhcnuggets-cwd /tools/mhcnuggets", command)
        self.assertIn("--netmhcpan-bin /tools/netMHCpan", command)

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

    def test_microbial_runtime_tmp_directory_is_created(self) -> None:
        tmp_dir = self.root / "nested" / "microbial-tmp"
        config = {"path": {"tmp_dir": str(tmp_dir)}}

        microbial._prepare_runtime_directories(config)

        self.assertTrue(tmp_dir.is_dir())
        self.assertEqual(config["path"]["tmp_dir"], str(tmp_dir.resolve()))

        with self.assertRaisesRegex(ValueError, "tmp_dir must be set"):
            microbial._prepare_runtime_directories({"path": {"tmp_dir": ""}})

    def test_cryptic_default_and_native_dispatch(self) -> None:
        config = self.cryptic_config()
        paths = self.cryptic_paths()
        tool = MagicMock()
        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T", config, paths, tool)
        self.assertEqual(run_cmd.call_count, 3)
        self.assertTrue(run_cmd.call_args_list[0].args[2][1].endswith("08-orf_genome_annotation.py"))
        self.assertTrue(any(str(value).endswith("/07-orf_genome_annotation") for value in run_cmd.call_args_list[0].args[2]))
        self.assertTrue(run_cmd.call_args_list[1].args[2][1].endswith("08-orf_filter.py"))
        self.assertTrue(any(str(value).endswith("/08-orf_filter") for value in run_cmd.call_args_list[1].args[2]))
        legacy_command = run_cmd.call_args.args[2]
        self.assertTrue(legacy_command[1].endswith("07-hla_binding_pred.py"))
        self.assertTrue(any(str(value).endswith("/09-hla_binding_pred") for value in legacy_command))
        self.assertTrue(any(str(value).endswith("CRYPTIC-T.aeSEPs.orf_filtered.pep") for value in legacy_command))

        config["others"].update(
            {
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_step_name": "09-hla_binding_pred_mimicneoai_test",
                "binding_prediction_preset": "fast",
                "binding_prediction_max_task_rows": 1234,
                "binding_prediction_force_large_samples": True,
            }
        )
        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T", config, paths, tool)
        native_command = run_cmd.call_args.args[2]
        self.assertTrue(native_command[1].endswith("07-hla_binding_pred_mimicneoai.py"))
        self.assertTrue(any(str(value).endswith("/09-hla_binding_pred_mimicneoai_test") for value in native_command))
        self.assertTrue(any(str(value).endswith("CRYPTIC-T.aeSEPs.orf_filtered.pep") for value in native_command))
        self.assertIn("--max-task-rows", native_command)
        self.assertIn("1234", native_command)
        self.assertIn("--preset", native_command)
        self.assertIn("fast", native_command)
        self.assertIn("--force-large-samples", native_command)
        self.assertIn("--netmhcpan-bin", native_command)
        self.assertIn("/tools/netMHCpan", native_command)

        config = self.cryptic_config()
        config["others"].update(
            {
                "cryptic_core_qc": True,
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_step_name": "09-hla_binding_pred_mimicneoai_core",
                "binding_prediction_preset": "fast",
                "allow_missing_human_reference": True,
            }
        )
        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T,CRYPTIC-N", config, paths, tool)
        self.assertEqual(run_cmd.call_count, 4)
        core_command = run_cmd.call_args_list[2].args[2]
        self.assertTrue(core_command[1].endswith("cryptic_core_qc.py"))
        self.assertIn("--matched-control-sample", core_command)
        self.assertIn("CRYPTIC-N", core_command)
        self.assertIn("--allow-missing-human-reference", core_command)
        binding_command = run_cmd.call_args_list[3].args[2]
        self.assertTrue(binding_command[1].endswith("07-hla_binding_pred_mimicneoai.py"))
        self.assertIn("--input-mode", binding_command)
        self.assertIn("peptide-core", binding_command)
        self.assertTrue(any(str(value).endswith("/08b.CrypticCoreQC_v1.0/cryptic_peptide_core.fasta") for value in binding_command))
        self.assertFalse(any(str(value).endswith("CRYPTIC-T.aeSEPs.orf_filtered.pep") for value in binding_command))

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
        paths = {
            "path": {
                "common": {
                    "BINDING_PREDICTORS": {
                        "NETMHCPAN_BIN": "/tools/netMHCpan"
                    }
                }
            }
        }
        tool = MagicMock()
        with patch.object(microbial_peptides, "pvacbind") as legacy:
            microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, paths, tool)
        legacy.assert_not_called()
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("hla_binding_pred_mimicneoai.py", command)
        self.assertIn("08.MicrobialPeptidesBindingPrediction_mimicneoai", command)
        self.assertIn("--preset fast", command)
        self.assertIn("--mhc-i-lengths 8,9,10,11", command)
        self.assertIn("--mhc-ii-lengths 13,14,15,16,17", command)
        self.assertNotIn("NNalign", command)
        self.assertIn("--netmhcpan-bin /tools/netMHCpan", command)

        config["others"]["binding_prediction_backend"] = "pvactools"
        tool.reset_mock()
        with patch.object(microbial_peptides, "pvacbind") as legacy:
            microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, paths, tool)
        legacy.assert_called_once()
        tool.exec_cmd.assert_not_called()

        config["others"].update(
            {
                "binding_prediction_backend": "mimicneoai",
                "binding_prediction_step_name": "08.MicrobialPeptidesBindingPrediction_mimicneoai_test",
                "binding_prediction_preset": "fast",
                "binding_prediction_max_task_rows": 4321,
                "binding_prediction_force_large_samples": True,
            }
        )
        tool.reset_mock()
        microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, paths, tool)
        command = tool.exec_cmd.call_args.args[0]
        self.assertIn("hla_binding_pred_mimicneoai.py", command)
        self.assertIn("08.MicrobialPeptidesBindingPrediction_mimicneoai_test", command)
        self.assertIn("--max-task-rows 4321", command)
        self.assertIn("--preset fast", command)
        self.assertIn("--force-large-samples", command)
        self.assertIn("--netmhcpan-bin /tools/netMHCpan", command)

        config["others"]["binding_prediction_step_name"] = "../invalid"
        with self.assertRaisesRegex(ValueError, "single directory name"):
            microbial_peptides.MicrobialPeptidesBindingPrediction(sample, config, paths, tool)

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
                "common": {
                    "PVACTOOLS": "/tools/pvactools.sif",
                    "BINDING_PREDICTORS": {
                        "NETMHCPAN_BIN": "/tools/netMHCpan"
                    },
                },
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
