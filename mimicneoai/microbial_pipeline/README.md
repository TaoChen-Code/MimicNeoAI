# Microbial Antigen Pipeline

Identify microbial peptides from host sequencing data and run binding/immunogenicity scoring.

## Overview
1. QC
2. Host depletion (hg38 + T2T)
3. Vector decontamination
4. Microbial taxa quantification
5. Microbial peptide identification
6. HLA typing
7. Binding prediction

## Installation

```bash
pip install mimicneoai
# or
conda install -c conda-forge -c bioconda mimicneoai
```

## External Dependencies (in `PATH`)

| Tool | Recommended/validated version | Purpose |
|---|---|---|
| `fastp` | `0.22.0` | FASTQ QC |
| `bwa` | `0.7.17` | Host/vector alignment |
| `samtools` | `1.5` | BAM/SAM/FASTQ processing |
| `java` | `17` | Run GATK PathSeq JAR |
| `blastx` | `2.15.0+` | Microbial peptide identification |
| `bowtie2` | `2.4.1` | HLA typing pre-alignment |
| `hlahd.sh` | `1.7.0` | HLA typing |
| `apptainer` | `1.4.2` | Run pVACtools container |
| `pVACtools` (container) | `4.2.1` | Binding prediction (`pvacbind`) |
| `GATK` (jar) | `4.6.0.0` | PathSeq |
| `awk` | system tool | FASTQ header normalization in HLA typing |

## Database and Paths

Shared documentation:
- [`mimicneoai/configures/Database_and_Paths.md`](../configures/Database_and_Paths.md)

## Configuration

Use the canonical template directly:
- [`mimicneoai/configures/microbial_configure.yaml`](../configures/microbial_configure.yaml)

Important: field names in this template are the source of truth. In particular, use current toggles such as:
- `others.run_host_depletion`
- `others.run_vector_decontamination`
- `others.run_pathseq`
- `others.run_microbial_peptide_identification`
- `others.run_hla_typing`
- `others.run_binding_prediction`

Matched-normal microbial mode uses the same sample syntax as the mutation
pipeline:

```yaml
others:
  tumor_with_matched_normal: true

samples:
  - Tumor_1,Normal_1
```

Each `Tumor,Normal` pair is parsed once as a task unit. Both samples run through
protein-hit identification, HLA typing is run only for the tumor, and binding is
run only on the tumor peptide Core after exact matched-normal peptide
subtraction. Single-sample mode remains the default when
`tumor_with_matched_normal` is false.

## Run

```bash
# Method 1
python -m mimicneoai.microbial_pipeline.microbial -c /path/to/microbial_configure.yaml

# Method 2 (unified CLI)
mimicneoai microbial -c /path/to/microbial_configure.yaml
```

## Output Structure

Pipeline outputs are written under:

```text
<output_dir>/Microbial/<sample>/
├── 00.QC
├── 01.HostSequencesRemovingStep1
├── 02.HostSequencesRemovingStep2
├── 03.VectorContaminationRemoving
├── 04.MicrobialTaxaQuantificationStep1
├── 05.MicrobialTaxaQuantificationStep2
├── 06.MicrobialPeptidesIdentification
├── 06b.MicrobialProteinCoreQC_v1.0
├── 07.HlaTyping
├── 08.MicrobialPeptidesBindingPrediction_mimicneoai
└── 09.ImmunogenicityPrediction_mimicneoai
```

`06.MicrobialPeptidesIdentification` writes both legacy and normalized
protein-hit products:

- `<sample>.blastx.filtered` and `<sample>.peptide.fasta` are retained for
  backward compatibility.
- `<sample>.protein_hits.filtered.tsv` is the normalized BLASTX/DIAMOND
  protein-hit table for newer downstream microbial Core construction.
- `<sample>.protein_hits.excluded.tsv` records fail-closed exclusions such as
  missing/failed coverage, non-100% identity, catalog mismatch, and
  noncanonical parent sequences.
- `<sample>.protein_hits.qc_summary.tsv` records stagewise protein-hit counts.

In matched-normal mode, `06b.MicrobialProteinCoreQC_v1.0` writes the paired Core:

- `microbial_parent_core.tsv` and `microbial_parent_excluded.tsv` retain parent
  protein-hit traceability.
- `microbial_peptide_core.tsv` is the tumor-only, matched-normal-depleted
  peptide space.
- `microbial_peptide_core_hla_i.fasta`,
  `microbial_peptide_core_hla_ii.fasta`, and `microbial_peptide_core.fasta`
  are already tiled candidate peptides. Binding consumes
  `microbial_peptide_core.fasta` in peptide-Core mode and does not tile it
  again.
- `matched_normal_peptide.tsv`, `microbial_peptide_parent_map.tsv`,
  `stagewise_qc.tsv`, and `run_manifest.json` provide subtraction and
  provenance QC.

## Notes

- The pipeline is resumable; existing non-empty outputs are skipped.
- If a step fails mid-way, delete the incomplete step directory and rerun.
- `others.binding_prediction_backend` defaults to `mimicneoai` with
  `others.binding_prediction_preset: fast`. This estimates task scale before
  materializing the task table and then uses EL-rank Stage 1 routing before
  formal local binding prediction; see the
  [native binding backend documentation](../functions/binding_prediction/README.md).
- Set `others.binding_prediction_preset: full` for one-stage multialgorithm
  prediction, or set `others.binding_prediction_backend: pvactools` to run the
  legacy pVACbind workflow.
- Set `others.run_immunogenicity_prediction: true` to score peptide-HLA rows
  after binding. Oversized samples that are skipped by the binding scale gate
  are also skipped for immunogenicity scoring and are not labeled negative.
- Immunogenicity inference requires a Python environment with PyTorch and
  scikit-learn. Use `others.immunogenicity_python_bin` or the
  `MIMICNEOAI_IMMUNOGENICITY_PYTHON_BIN` environment variable to select a CPU
  or GPU runtime explicitly; otherwise the current pipeline Python is used.
