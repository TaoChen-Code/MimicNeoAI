from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from mimicneoai.cryptic_pipeline.scripts.cryptic_core_qc import build_core, load_human_reference_matches
from mimicneoai.cryptic_pipeline.scripts.cryptic_junction_qc import manifest_critical_parameter_digest


class CrypticCoreQCTest(unittest.TestCase):
    codons = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "TTA",
        "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
        "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
    }

    def _nt_for_peptide(self, peptide: str) -> str:
        return "".join(self.codons[aa] for aa in peptide)

    def _reverse_complement(self, seq: str) -> str:
        return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def _sha256_file(self, path: Path) -> str:
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    def _write_star_pair(self, root: Path, tumor_sj: Path, normal_sj: Path) -> Path:
        tumor_sha = self._sha256_file(tumor_sj)
        normal_sha = self._sha256_file(normal_sj)
        star_cmd = (
            "STAR --runThreadN 8 --genomeDir {index} --readFilesIn R1.fq.gz R2.fq.gz "
            "--readFilesCommand zcat --outSAMtype BAM Unsorted "
            "--outFilterMultimapNmax 20 --alignSJoverhangMin 8 "
            "--alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 "
            "--alignMatesGapMax 1000000 --outFileNamePrefix out/"
        ).format(index=root / "STAR-index")
        tumor_log = root / "tumor.Log.out"
        normal_log = root / "normal.Log.out"
        tumor_log.write_text(f"##### Command Line:\n{star_cmd}\n")
        normal_log.write_text(f"##### Command Line:\n{star_cmd}\n")
        tumor_log_sha = self._sha256_file(tumor_log)
        normal_log_sha = self._sha256_file(normal_log)
        critical_parameters = {
            "--readFilesCommand": ["zcat"],
            "--outSAMtype": ["BAM", "Unsorted"],
            "--outFilterMultimapNmax": ["20"],
            "--alignSJoverhangMin": ["8"],
            "--alignSJDBoverhangMin": ["1"],
            "--alignIntronMin": ["20"],
            "--alignIntronMax": ["1000000"],
            "--alignMatesGapMax": ["1000000"],
        }
        critical_digest = manifest_critical_parameter_digest(critical_parameters)
        tumor_manifest = root / "tumor_star.json"
        normal_manifest = root / "normal_star.json"
        source_manifest = root / "normal_source_star_alignment.manifest.json"
        source_manifest.write_text(json.dumps({
            "status": "completed",
            "input_signature": {
                "alignment_role": "control",
                "sample": "CTRL",
                "tumor_sample": "TEST",
                "control_sample": "CTRL",
            },
        }))
        source_manifest_sha = self._sha256_file(source_manifest)
        common = {
            "policy_version": "star_provenance_freeze_v1.0",
            "tumor_sample": "TEST",
            "normal_sample": "CTRL",
            "star_version": "STAR_2.5.3a_modified",
            "star_index": str(root / "STAR-index"),
        }
        tumor_manifest.write_text(json.dumps({
            **common,
            "status": "frozen_legacy_complete",
            "critical_parameter_compatibility": "compatible",
            "critical_star_parameters": critical_parameters,
            "critical_star_parameter_digest": critical_digest,
            "current_outputs": {
                "SJ.out.tab": {
                    "path": str(tumor_sj),
                    "sha256": tumor_sha,
                    "size_bytes": tumor_sj.stat().st_size,
                },
                "Log.out": {
                    "path": str(tumor_log),
                    "sha256": tumor_log_sha,
                    "size_bytes": tumor_log.stat().st_size,
                },
            },
        }))
        normal_manifest.write_text(json.dumps({
            **common,
            "status": "relocated_complete",
            "critical_parameter_compatibility": "compatible",
            "critical_star_parameters": critical_parameters,
            "critical_star_parameter_digest": critical_digest,
            "source_manifest_path": str(source_manifest),
            "source_manifest_sha256": source_manifest_sha,
            "current_outputs": {
                "SJ.out.tab": {
                    "path": str(normal_sj),
                    "sha256": normal_sha,
                    "size_bytes": normal_sj.stat().st_size,
                },
                "Log.out": {
                    "path": str(normal_log),
                    "sha256": normal_log_sha,
                    "size_bytes": normal_log.stat().st_size,
                },
            },
        }))
        star_pairs = root / "star_pairs.tsv"
        star_pairs.write_text(
            "tumor_sample\tnormal_sample\ttumor_sj_path\tnormal_sj_path\t"
            "tumor_sj_sha256\tnormal_sj_sha256\ttumor_star_manifest\tnormal_star_manifest\t"
            "critical_parameter_compatibility\ttumor_critical_parameter_digest\tnormal_critical_parameter_digest\n"
            f"TEST\tCTRL\t{tumor_sj}\t{normal_sj}\t{tumor_sha}\t{normal_sha}\t"
            f"{tumor_manifest}\t{normal_manifest}\tcompatible\t{critical_digest}\t{critical_digest}\n"
        )
        return star_pairs

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

    def _write_coordinate_inputs(self, root: Path, parent_sequences: dict[str, str]) -> dict[str, Path]:
        import pysam

        bed12 = root / "orf.noUnmap.noSup.bed12"
        bam = root / "orf2genome.noUnmap.noSup.bam"
        cds = root / "sample.SEPs.cds.fa"
        ref = root / "GRCh38.test.fa"
        contig = ["N"] * 20000

        def place(start0: int, seq: str) -> None:
            contig[start0 : start0 + len(seq)] = list(seq)

        with bed12.open("w") as handle:
            for idx, (parent_id, peptide) in enumerate(parent_sequences.items()):
                start = 100 + idx * 1000
                cds_seq = self._nt_for_peptide(peptide)
                size = len(cds_seq)
                strand = "-"
                block_sizes = [size]
                block_starts = [0]
                if parent_id == "ENST_NC.1.p1":
                    strand = "+"
                    block_sizes = [18, size - 18]
                    block_starts = [0, 100]
                    place(start, cds_seq[:18])
                    place(start + 100, cds_seq[18:])
                    chrom_end = start + 100 + block_sizes[-1]
                elif parent_id.startswith("TRINITY"):
                    place(start, self._reverse_complement(cds_seq))
                    chrom_end = start + size
                else:
                    strand = "+"
                    place(start, cds_seq)
                    chrom_end = start + size
                handle.write(
                    f"chr1\t{start}\t{chrom_end}\t{parent_id}\t30\t{strand}\t"
                    f"{start}\t{chrom_end}\t0\t{len(block_sizes)}\t"
                    f"{','.join(str(v) for v in block_sizes)},\t"
                    f"{','.join(str(v) for v in block_starts)},\n"
                )

        ref.write_text(">chr1\n" + "".join(contig) + "\n")
        pysam.faidx(str(ref))

        with cds.open("w") as handle:
            for parent_id, peptide in parent_sequences.items():
                handle.write(f">{parent_id}\n{self._nt_for_peptide(peptide)}\n")

        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000000}]}
        with pysam.AlignmentFile(bam, "wb", header=header) as out:
            for idx, (parent_id, peptide) in enumerate(parent_sequences.items()):
                start = 100 + idx * 1000
                size = len(peptide) * 3
                rec = pysam.AlignedSegment()
                rec.query_name = parent_id
                rec.query_sequence = "A" * size
                rec.flag = 16 if parent_id.startswith("TRINITY") else 0
                rec.reference_id = 0
                rec.reference_start = start
                rec.mapping_quality = 30
                if parent_id == "ENST_NC.1.p1":
                    rec.cigar = [(0, 18), (3, 82), (0, size - 18)]
                else:
                    rec.cigar = [(0, size)]
                out.write(rec)
                if parent_id.endswith("DUP.1.p1"):
                    sec = pysam.AlignedSegment()
                    sec.query_name = parent_id
                    sec.query_sequence = "A" * size
                    sec.flag = 256
                    sec.reference_id = 0
                    sec.reference_start = start + 200
                    sec.mapping_quality = 30
                    sec.cigar = [(0, size)]
                    out.write(sec)
                    sup = pysam.AlignedSegment()
                    sup.query_name = parent_id
                    sup.query_sequence = "A" * size
                    sup.flag = 2048
                    sup.reference_id = 0
                    sup.reference_start = start + 400
                    sup.mapping_quality = 30
                    sup.cigar = [(0, size)]
                    out.write(sup)
        return {"bed12": bed12, "bam": bam, "cds": cds, "ref": ref}

    def _args(self, paths: dict[str, Path], outdir: Path) -> argparse.Namespace:
        return argparse.Namespace(
            sample="TEST",
            policy_version="cryptic_core_qc_v1.0",
            ae_seps_fasta=str(paths["fasta"]),
            aeseps_annotation=str(paths["annot"]),
            orf_filtered_fasta=str(paths["fasta"]),
            orf_final=str(paths["orf_final"]),
            orf_bed12="",
            orf_bam="",
            orf_cds_fasta="",
            outdir=str(outdir),
            human_proteome_fasta=str(paths["human"]),
            allow_missing_human_reference=False,
            matched_control_sample="CTRL",
            reference_genome_fasta="",
            reference_gtf="",
            reference_lnc_gtf="",
            reference_build="GRCh38",
            strandedness="reverse",
            min_tpm_tumor=5.0,
            max_tpm_ctrl=0.5,
            min_log2fc=4.0,
            mhc_i_lengths="8",
            mhc_ii_lengths="13",
            candidate_selection_mode="all",
            max_hla_i_peptides=None,
            max_hla_ii_peptides=None,
            coordinate_min_mapq=20,
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

    def test_ranked_cap_limits_unique_peptides_and_records_deferred_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            outdir = root / "ranked"
            args = self._args(paths, outdir)
            args.candidate_selection_mode = "ranked_cap"
            args.max_hla_i_peptides = 1
            args.max_hla_ii_peptides = 1

            manifest = build_core(args)

            self.assertEqual(manifest["candidate_selection"]["mode"], "ranked_cap")
            self.assertEqual(manifest["stage_counts"]["hla_i_unique_peptide_core_records"], 1)
            self.assertEqual(manifest["stage_counts"]["hla_ii_unique_peptide_core_records"], 1)
            self.assertTrue((outdir / "cryptic_parent_ranked.tsv").exists())
            self.assertTrue((outdir / "cryptic_peptide_selection_evidence.tsv").exists())
            self.assertTrue((outdir / "cryptic_deferred_parent.tsv").exists())
            self.assertTrue((outdir / "cryptic_peptide_deferred.tsv").exists())
            evidence = pd.read_csv(outdir / "cryptic_peptide_selection_evidence.tsv", sep="\t")
            self.assertIn("human_reference_peptide_match", set(evidence["selection_reason"].dropna().astype(str)))
            self.assertIn("selected_for_binding", set(evidence["selection_status"]))
            parent_map = pd.read_csv(outdir / "cryptic_peptide_parent_map.tsv", sep="\t")
            selected = parent_map[parent_map["peptide_record_id"].fillna("").astype(str).ne("")]
            self.assertEqual(selected[["mhc_class", "peptide"]].drop_duplicates().shape[0], 2)
            stage = pd.read_csv(outdir / "stagewise_qc.tsv", sep="\t").set_index("stage")["count"].to_dict()
            self.assertGreater(stage["peptide_deferred_window_rows"], 0)
            self.assertEqual(
                int(stage["peptide_excluded_window_rows"]),
                int(parent_map["human_reference_peptide_status"].eq("human_reference_peptide_match").sum()),
            )
            deferred = pd.read_csv(outdir / "cryptic_peptide_deferred.tsv", sep="\t")
            self.assertTrue(deferred["peptide_core_status"].eq("deferred").all())
            self.assertNotIn("not_selected_due_to_analysis_cap", set(parent_map["peptide_qc_reasons"].dropna().astype(str)))

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

    def test_v11_writes_coordinate_sidecars_and_marks_secondary_ambiguity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            coordinate_inputs = self._write_coordinate_inputs(
                root,
                {
                    "ENST_NC.1.p1": "ACDEFGHIKLMNP",
                    "TRINITY_X|NOVEL|.|..p1": "KLMNPQRSTVWY",
                    "ENST_DUP.1.p1": "ACDEFGHIKLMNP",
                },
            )
            outdir = root / "core_v11"
            args = self._args(paths, outdir)
            args.policy_version = "cryptic_core_qc_v1.1"
            args.orf_bed12 = str(coordinate_inputs["bed12"])
            args.orf_bam = str(coordinate_inputs["bam"])
            args.orf_cds_fasta = str(coordinate_inputs["cds"])
            args.reference_genome_fasta = str(coordinate_inputs["ref"])

            manifest = build_core(args)

            self.assertEqual(manifest["policy_version"], "cryptic_core_qc_v1.1")
            self.assertIn("cryptic_parent_coordinates.tsv", manifest["output_signature"])
            self.assertIn("cryptic_parent_orfcds.tsv", manifest["output_signature"])
            self.assertIn("cryptic_peptide_genomic_footprint.tsv", manifest["output_signature"])
            parent_coords = pd.read_csv(outdir / "cryptic_parent_coordinates.tsv", sep="\t")
            status_by_parent = parent_coords.set_index("parent_record_id")["coordinate_mapping_status"].to_dict()
            self.assertEqual(status_by_parent["ENST_NC.1.p1"], "coordinate_evaluable")
            self.assertEqual(status_by_parent["TRINITY_X|NOVEL|.|..p1"], "coordinate_evaluable")
            self.assertEqual(status_by_parent["ENST_DUP.1.p1"], "ambiguous_candidate_mapping")
            self.assertEqual(
                parent_coords.set_index("parent_record_id").loc[
                    "ENST_NC.1.p1", "reference_genomic_translation_status"
                ],
                "pass",
            )
            self.assertEqual(
                parent_coords.set_index("parent_record_id").loc[
                    "TRINITY_X|NOVEL|.|..p1", "reference_genomic_translation_status"
                ],
                "pass",
            )
            self.assertEqual(
                int(parent_coords.set_index("parent_record_id").loc[
                    "ENST_DUP.1.p1", "supplementary_alignment_count"
                ]),
                1,
            )
            footprints = pd.read_csv(outdir / "cryptic_peptide_genomic_footprint.tsv", sep="\t")
            self.assertIn("coordinate_evaluable", set(footprints["candidate_coordinate_status"]))
            self.assertIn("ambiguous_candidate_mapping", set(footprints["candidate_coordinate_status"]))
            self.assertEqual(manifest["stage_counts"]["candidate_parents_coordinate_evaluable"], 2)
            self.assertEqual(manifest["stage_counts"]["candidate_parents_coordinate_not_evaluable"], 1)
            self.assertEqual(manifest["stage_counts"]["candidate_parents_reference_translation_pass"], 3)

    def test_v11_junction_qc_filters_before_peptide_core(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            coordinate_inputs = self._write_coordinate_inputs(
                root,
                {
                    "ENST_NC.1.p1": "ACDEFGHIKLMNP",
                    "TRINITY_X|NOVEL|.|..p1": "KLMNPQRSTVWY",
                    "ENST_DUP.1.p1": "ACDEFGHIKLMNP",
                },
            )
            tumor_sj = root / "tumor.SJ.out.tab"
            normal_sj = root / "normal.SJ.out.tab"
            tumor_sj.write_text("chr1\t119\t200\t1\t1\t1\t2\t0\t20\n")
            normal_sj.write_text("chr1\t119\t200\t1\t1\t1\t1\t0\t20\n")
            star_pairs = self._write_star_pair(root, tumor_sj, normal_sj)
            outdir = root / "core_v11_junction"
            args = self._args(paths, outdir)
            args.policy_version = "cryptic_core_qc_v1.1"
            args.orf_bed12 = str(coordinate_inputs["bed12"])
            args.orf_bam = str(coordinate_inputs["bam"])
            args.orf_cds_fasta = str(coordinate_inputs["cds"])
            args.reference_genome_fasta = str(coordinate_inputs["ref"])
            args.junction_qc_enabled = True
            args.junction_policy_version = "junction_qc_v1.0"
            args.star_pair_inputs = str(star_pairs)
            args.primary_min_tumor_unique_reads = 2
            args.junction_sensitivity_thresholds = "1,2,3,5"

            manifest = build_core(args)

            core_parents = pd.read_csv(outdir / "cryptic_parent_core.tsv", sep="\t")
            excluded = pd.read_csv(outdir / "cryptic_parent_excluded.tsv", sep="\t")
            junction_summary = pd.read_csv(outdir / "cryptic_parent_junction_summary.tsv", sep="\t")
            self.assertEqual(set(core_parents["parent_record_id"]), {"ENST_NC.1.p1", "TRINITY_X|NOVEL|.|..p1"})
            self.assertIn("ENST_DUP.1.p1", set(excluded["parent_record_id"]))
            self.assertTrue((outdir / "cryptic_peptide_junction_evidence.tsv").exists())
            self.assertEqual(
                junction_summary.set_index("parent_record_id").loc["ENST_NC.1.p1", "junction_qc_status"],
                "all_required_junctions_tumor_unique_ge2",
            )
            self.assertEqual(manifest["junction_qc"]["enabled"], True)
            self.assertEqual(manifest["stage_counts"]["junction_parent_primary_core_eligible"], 2)
            self.assertIn("junction_qc_sha256", manifest["input_signature"]["junction_qc"])
            self.assertIn("pair_row_sha256", manifest["input_signature"]["junction_qc"]["star_pair_validation"])

            saved_path = outdir / "run_manifest.json"
            saved = json.loads(saved_path.read_text())
            saved["input_signature"]["junction_qc"]["junction_qc_sha256"] = "changed"
            saved_path.write_text(json.dumps(saved, indent=2))
            with self.assertRaises(ValueError):
                build_core(args)

    def test_v11_junction_qc_rejects_ranked_cap_until_refill_qc_is_implemented(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            coordinate_inputs = self._write_coordinate_inputs(root, {"ENST_NC.1.p1": "ACDEFGHIKLMNP"})
            tumor_sj = root / "tumor.SJ.out.tab"
            normal_sj = root / "normal.SJ.out.tab"
            tumor_sj.write_text("chr1\t119\t200\t1\t1\t1\t2\t0\t20\n")
            normal_sj.write_text("chr1\t119\t200\t1\t1\t1\t0\t0\t20\n")
            args = self._args(paths, root / "core_v11_ranked_junction")
            args.policy_version = "cryptic_core_qc_v1.1"
            args.orf_bed12 = str(coordinate_inputs["bed12"])
            args.orf_bam = str(coordinate_inputs["bam"])
            args.orf_cds_fasta = str(coordinate_inputs["cds"])
            args.reference_genome_fasta = str(coordinate_inputs["ref"])
            args.junction_qc_enabled = True
            args.junction_policy_version = "junction_qc_v1.0"
            args.star_pair_inputs = str(self._write_star_pair(root, tumor_sj, normal_sj))
            args.primary_min_tumor_unique_reads = 2
            args.junction_sensitivity_thresholds = "1,2,3,5"
            args.candidate_selection_mode = "ranked_cap"
            args.max_hla_i_peptides = 1
            args.max_hla_ii_peptides = 1
            with self.assertRaisesRegex(ValueError, "ranked_cap with junction_qc is disabled"):
                build_core(args)

    def test_v11_reference_translation_mismatch_is_not_coordinate_evaluable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._write_inputs(root)
            coordinate_inputs = self._write_coordinate_inputs(
                root,
                {
                    "ENST_NC.1.p1": "ACDEFGHIKLMNP",
                    "TRINITY_X|NOVEL|.|..p1": "KLMNPQRSTVWY",
                    "ENST_DUP.1.p1": "ACDEFGHIKLMNP",
                },
            )
            ref = Path(coordinate_inputs["ref"])
            text = ref.read_text()
            ref.write_text(text.replace("GCT", "AAA", 1))
            import pysam
            pysam.faidx(str(ref))

            outdir = root / "core_v11_mismatch"
            args = self._args(paths, outdir)
            args.policy_version = "cryptic_core_qc_v1.1"
            args.orf_bed12 = str(coordinate_inputs["bed12"])
            args.orf_bam = str(coordinate_inputs["bam"])
            args.orf_cds_fasta = str(coordinate_inputs["cds"])
            args.reference_genome_fasta = str(ref)

            manifest = build_core(args)
            parent_coords = pd.read_csv(outdir / "cryptic_parent_coordinates.tsv", sep="\t")
            row = parent_coords.set_index("parent_record_id").loc["ENST_NC.1.p1"]
            self.assertEqual(row["reference_genomic_translation_status"], "reference_translation_mismatch")
            self.assertIn("reference_translation_mismatch", row["coordinate_mapping_reasons"])
            self.assertEqual(row["rna_variant_coordinate_status"], "RNA_variant_aware_not_evaluated")
            self.assertEqual(manifest["stage_counts"]["candidate_parents_reference_translation_mismatch"], 1)


if __name__ == "__main__":
    unittest.main()
