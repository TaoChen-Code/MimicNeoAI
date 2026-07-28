from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from mimicneoai.cryptic_pipeline import cryptic


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ALIGNMENT_SCRIPT = PACKAGE_ROOT / "cryptic_pipeline" / "scripts" / "01-alignment.py"


def load_alignment_module():
    spec = importlib.util.spec_from_file_location("cryptic_star_alignment", ALIGNMENT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class MatchedNormalAlignmentDispatchTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.root = Path(self.tempdir.name)

    def _touch_raw_fastqs(self, *samples: str) -> None:
        for sample in samples:
            sample_dir = self.root / "input" / sample
            sample_dir.mkdir(parents=True, exist_ok=True)
            (sample_dir / f"{sample}.R1.fq.gz").write_text("r1\n", encoding="utf-8")
            (sample_dir / f"{sample}.R2.fq.gz").write_text("r2\n", encoding="utf-8")

    def _config(self, **others_updates) -> dict[str, object]:
        others = {
            "QC": True,
            "alignment": True,
            "alignment_control": False,
            "known": False,
            "novel": False,
            "salmon_quant": False,
            "salmon_quant_control": False,
            "hlatyping": False,
            "extract_aeseps": False,
            "orf_genome_annotation": False,
            "orf_filter": False,
            "cryptic_core_qc": False,
            "cryptic_external_normal_qc": False,
            "hla_binding_pred": False,
        }
        others.update(others_updates)
        return {
            "path": {"input_dir": str(self.root / "input"), "output_dir": str(self.root)},
            "args": {"threads": 2, "hla_binding_threads": 2},
            "others": others,
        }

    @staticmethod
    def _paths() -> dict[str, object]:
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

    def _commands(self, run_cmd) -> list[list[str]]:
        return [call.args[2] for call in run_cmd.call_args_list]

    def test_paired_alignment_control_dispatches_tumor_and_control_star(self) -> None:
        self._touch_raw_fastqs("CRYPTIC-T", "CRYPTIC-N")
        config = self._config(alignment_control=True)

        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T,CRYPTIC-N", config, self._paths(), MagicMock())

        commands = self._commands(run_cmd)
        star_commands = [cmd for cmd in commands if cmd[1].endswith("01-alignment.py")]
        self.assertEqual(len(star_commands), 2)
        self.assertIn("CRYPTIC-T", star_commands[0])
        self.assertIn("tumor", star_commands[0])
        self.assertIn("CRYPTIC-N", star_commands[1])
        self.assertIn("control", star_commands[1])
        self.assertTrue(any("--raw-fq1" in cmd for cmd in star_commands))
        self.assertEqual(run_cmd.call_args_list[2].kwargs["display_name"], "STAR alignment")
        self.assertEqual(run_cmd.call_args_list[3].kwargs["display_name"], "Control STAR alignment")

    def test_alignment_control_false_does_not_run_control_star(self) -> None:
        self._touch_raw_fastqs("CRYPTIC-T", "CRYPTIC-N")

        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T,CRYPTIC-N", self._config(), self._paths(), MagicMock())

        star_commands = [cmd for cmd in self._commands(run_cmd) if cmd[1].endswith("01-alignment.py")]
        self.assertEqual(len(star_commands), 1)
        self.assertIn("CRYPTIC-T", star_commands[0])
        self.assertNotIn("CRYPTIC-N", star_commands[0])

    def test_single_sample_never_runs_control_star(self) -> None:
        self._touch_raw_fastqs("CRYPTIC-T")
        config = self._config(alignment_control=False)

        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T", config, self._paths(), MagicMock())

        star_commands = [cmd for cmd in self._commands(run_cmd) if cmd[1].endswith("01-alignment.py")]
        self.assertEqual(len(star_commands), 1)
        self.assertIn("CRYPTIC-T", star_commands[0])

    def test_alignment_control_requires_existing_control_fastqs(self) -> None:
        self._touch_raw_fastqs("CRYPTIC-T")
        config = self._config(alignment_control=True)

        with (
            patch.object(cryptic, "_run_cmd") as run_cmd,
            self.assertRaisesRegex(FileNotFoundError, "control FASTQ"),
        ):
            cryptic._run_one_sample("CRYPTIC-T,CRYPTIC-N", config, self._paths(), MagicMock())
        run_cmd.assert_not_called()

    def test_control_star_does_not_make_normal_downstream_discovery(self) -> None:
        self._touch_raw_fastqs("CRYPTIC-T", "CRYPTIC-N")
        config = self._config(alignment_control=True, known=True, novel=True)

        with patch.object(cryptic, "_run_cmd") as run_cmd:
            cryptic._run_one_sample("CRYPTIC-T,CRYPTIC-N", config, self._paths(), MagicMock())

        downstream_commands = [
            cmd for cmd in self._commands(run_cmd) if cmd[1].endswith("02-lnc_sORF_pipeline.py")
        ]
        self.assertEqual(len(downstream_commands), 2)
        for cmd in downstream_commands:
            self.assertIn("CRYPTIC-T", cmd)
            self.assertNotIn("CRYPTIC-N", cmd)


class StarAlignmentCompletionTest(unittest.TestCase):
    def _signature_inputs(self, alignment, tmpdir: Path, sample: str) -> tuple[dict[str, object], list[str], Path, Path]:
        clean_r1 = tmpdir / f"{sample}.R1.QC.fq.gz"
        clean_r2 = tmpdir / f"{sample}.R2.QC.fq.gz"
        clean_r1.write_text("r1\n", encoding="utf-8")
        clean_r2.write_text("r2\n", encoding="utf-8")
        cmd = [
            "STAR",
            "--genomeDir", str(tmpdir / "star-index"),
            "--runThreadN", "2",
            "--outFileNamePrefix", str(tmpdir / "star" / sample),
            "--readFilesIn", str(clean_r1), str(clean_r2),
        ]
        signature = alignment.alignment_signature(
            sample=sample,
            alignment_role="control",
            tumor_sample="CRYPTIC-T",
            control_sample=sample,
            raw_r1="",
            raw_r2="",
            clean_r1=clean_r1,
            clean_r2=clean_r2,
            genome_dir=str(tmpdir / "star-index"),
            star_bin="STAR",
            threads="2",
            command=cmd,
        )
        return signature, cmd, clean_r1, clean_r2

    def _write_complete_outputs(self, star_dir: Path, sample: str) -> None:
        for suffix in ("Aligned.out.bam", "SJ.out.tab", "Log.final.out", "Log.out"):
            (star_dir / f"{sample}{suffix}").write_text(f"{suffix}\n", encoding="utf-8")

    def test_resume_requires_bam_sj_final_log_and_log_out(self) -> None:
        alignment = load_alignment_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            star_dir = Path(tmpdir)
            sample = "CRYPTIC-N"
            (star_dir / f"{sample}Aligned.out.bam").write_text("bam\n", encoding="utf-8")
            self.assertFalse(alignment.star_outputs_complete(star_dir, sample))

            (star_dir / f"{sample}SJ.out.tab").write_text("sj\n", encoding="utf-8")
            self.assertFalse(alignment.star_outputs_complete(star_dir, sample))

            (star_dir / f"{sample}Log.final.out").write_text("log\n", encoding="utf-8")
            self.assertFalse(alignment.star_outputs_complete(star_dir, sample))

            (star_dir / f"{sample}Log.out").write_text("log\n", encoding="utf-8")
            self.assertTrue(alignment.star_outputs_complete(star_dir, sample))

    def test_complete_output_requires_manifest_signature_match(self) -> None:
        alignment = load_alignment_module()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            sample = "CRYPTIC-N"
            star_dir = tmpdir / "star"
            star_dir.mkdir()
            manifest = star_dir / "star_alignment.manifest.json"
            self._write_complete_outputs(star_dir, sample)
            signature, cmd, clean_r1, clean_r2 = self._signature_inputs(alignment, tmpdir, sample)

            with self.assertRaisesRegex(RuntimeError, "manifest is missing"):
                alignment.existing_star_state(star_dir, sample, manifest, signature)

            alignment.write_manifest(
                manifest,
                sample=sample,
                alignment_role="control",
                tumor_sample="CRYPTIC-T",
                control_sample=sample,
                raw_r1="",
                raw_r2="",
                clean_r1=clean_r1,
                clean_r2=clean_r2,
                genome_dir=str(tmpdir / "star-index"),
                star_bin="STAR",
                threads="2",
                star_dir=star_dir,
                command=cmd,
                status="completed",
            )
            self.assertEqual(
                alignment.existing_star_state(star_dir, sample, manifest, signature),
                "complete",
            )

            mismatch = dict(signature)
            mismatch["threads"] = 4
            with self.assertRaisesRegex(RuntimeError, "signature does not match"):
                alignment.existing_star_state(star_dir, sample, manifest, mismatch)

    def test_partial_output_fails_closed_instead_of_overwriting(self) -> None:
        alignment = load_alignment_module()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            sample = "CRYPTIC-N"
            star_dir = tmpdir / "star"
            star_dir.mkdir()
            (star_dir / f"{sample}Aligned.out.bam").write_text("bam\n", encoding="utf-8")
            signature, _, _, _ = self._signature_inputs(alignment, tmpdir, sample)

            with self.assertRaisesRegex(RuntimeError, "Incomplete STAR output"):
                alignment.existing_star_state(
                    star_dir,
                    sample,
                    star_dir / "star_alignment.manifest.json",
                    signature,
                )


if __name__ == "__main__":
    unittest.main()
