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
├── 07.HlaTyping
└── 08.MicrobialPeptidesBindingPrediction
```

## Notes

- The pipeline is resumable; existing non-empty outputs are skipped.
- If a step fails mid-way, delete the incomplete step directory and rerun.
- `others.binding_prediction_backend` defaults to `pvactools`. The optional
  `mimicneoai` backend estimates task scale before materializing the task table;
  see the [native binding backend documentation](../functions/binding_prediction/README.md).
