# Cryptic (sORF-Encoded) Antigen Pipeline

Discover and prioritize sORF-encoded peptides from known and novel transcripts.

## Overview
1. QC
2. STAR alignment
3. Known/novel lncRNA-sORF discovery
4. Tumor/control quantification
5. HLA typing
6. Aberrantly expressed sORF peptide extraction
7. ORF genome annotation and ORF-level filtering
8. Cryptic Core QC
9. External-normal exact sequence QC
10. Binding prediction

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
| `STAR` | `2.5.3a` | RNA alignment |
| `samtools` | `1.5` | BAM/SAM/FASTQ processing |
| `stringtie` | `3.0.1` | Transcript assembly (novel branch) |
| `gffcompare` | `0.12.10` | Transcript annotation comparison |
| `gffread` | `0.12.7` | GTF/FASTA conversion |
| `minimap2` | `2.30-r1287` | Contig-to-reference alignment |
| `bcftools` | `1.11` | Variant calling/filtering/consensus (known branch) |
| `TransDecoder.LongOrfs` | `5.5.0` | ORF calling |
| `salmon` | `1.10.0` | Expression quantification |
| `bowtie2` | `2.4.1` | HLA typing pre-alignment |
| `hlahd.sh` | `1.7.0` | HLA typing |
| `apptainer` | `1.4.2` | Trinity/pVACtools container execution |
| `pVACtools` (container) | `4.2.1` | Binding prediction (`pvacbind`) |
| `awk` | system tool | FASTQ header normalization in HLA typing |

Optional (only when `others.trinity_mode: native`):
- `Trinity`

## Database and Paths

Shared documentation:
- [`mimicneoai/configures/Database_and_Paths.md`](../configures/Database_and_Paths.md)

## Configuration

Use the canonical template directly:
- [`mimicneoai/configures/cryptic_configure.yaml`](../configures/cryptic_configure.yaml)

Key switches are defined in this template (for example `others.known`, `others.novel`, `others.salmon_quant_control`, `others.hla_binding_pred`). Please follow template key names exactly.

## Run

```bash
# Method 1
python -m mimicneoai.cryptic_pipeline.cryptic -c /path/to/cryptic_configure.yaml

# Method 2 (unified CLI)
mimicneoai cryptic -c /path/to/cryptic_configure.yaml
```

## Output Structure

Pipeline outputs are written under:

```text
<output_dir>/Cryptic/<tumor_sample>/
├── 00-clean
├── 01-star
├── 02-known
├── 03-novel
├── 04-salmon_quant
├── 05-hla_typing
├── 06-aeSEPs
├── 07-orf_genome_annotation
├── 08-orf_filter
├── 08b-cryptic_core_qc
├── 08c-external_normal_qc
├── 09-hla_binding_pred_mimicneoai
├── 10-immunogenicity_prediction_mimicneoai
└── 023-shared
```

Notable subfolders:
- `04-salmon_quant/salmon_index`, `salmon_quant`, `salmon_quant_control`
- `07-orf_genome_annotation`: maps selected ORF/CDS records back to the reference genome
- `08-orf_filter`: writes the ORF-filtered aeSEP FASTA used by binding prediction
- `08b-cryptic_core_qc`: materializes the strict pre-binding Cryptic Core
  with parent sidecars, HLA-I/HLA-II peptide-core FASTA files, stagewise counts,
  and an input/config manifest.
- `08c-external_normal_qc`: applies frozen external-normal sequence
  resources using exact peptide sequence plus HLA class. It keeps the complete
  source-Core landscape and writes the tumor-restricted primary Core FASTA used
  by downstream binding.
- `09-hla_binding_pred/<tumor_sample>/pvacbind` for the pVACbind backend, or
  `09-hla_binding_pred_mimicneoai/<tumor_sample>/` for the MimicNeoAI backend

## Notes

- Tumor/control samples should be provided as `Tumor,Control` in `samples`.
- The pipeline is resumable when manifest, input, code, configuration, and
  output signatures match. Empty Core FASTA files are valid for zero-candidate
  samples and can be reused through the manifest.
- `others.orf_genome_annotation` and `others.orf_filter` default to enabled;
  binding prediction uses the ORF-filtered aeSEP FASTA when `orf_filter` is enabled.
- `others.cryptic_core_qc` defaults to enabled in the template. It accepts only
  `novel` and `noncoding` aeSEP sources, requires ORF-filtered parent records,
  preserves excluded rows, and writes pre-tiled peptide-core FASTA files for
  HLA-I 8-11 aa and HLA-II 13-17 aa.
- `others.cryptic_external_normal_qc` defaults to disabled in the generic
  template because the frozen resource package is project-specific. Formal
  project configs should enable it and provide `external_normal_resources`, or
  rely on populated `paths.yaml` resource entries. Missing or hash-mismatched
  resources fail closed. `allow_missing_external_normal_resources: true` is
  exploratory only and cannot be routed into binding.
- `others.binding_prediction_backend` defaults to `mimicneoai` in the template.
  The local backend estimates task scale before materializing the task table;
  see the [native binding backend documentation](../functions/binding_prediction/README.md).
- For legacy pVACtools reruns, disable `others.cryptic_core_qc` or explicitly
  provide a parent FASTA. The strict peptide-core output is routed only to the
  MimicNeoAI backend to avoid accidental second-pass peptide tiling.
- When `08c` is enabled, binding consumes
  `08c-external_normal_qc/cryptic_tumor_restricted_primary_core.fasta`.
  It does not fall back to `08b`, `06-aeSEPs`, or the ORF-filtered parent FASTA.
- With `others.binding_prediction_backend: mimicneoai`, set
  `others.binding_prediction_preset: full` for one-stage multialgorithm
  prediction or `fast` for EL-rank Stage 1 routing before formal local binding
  prediction. Omitting the preset preserves explicit length and algorithm
  settings in the YAML.
- Set `others.run_immunogenicity_prediction: true` to score peptide-HLA rows
  after binding. If binding is skipped by the scale gate, immunogenicity scoring
  is skipped and recorded in the summary rather than being treated as negative.
- Immunogenicity inference requires a Python environment with PyTorch and
  scikit-learn. Use `others.immunogenicity_python_bin` or the
  `MIMICNEOAI_IMMUNOGENICITY_PYTHON_BIN` environment variable to select a CPU
  or GPU runtime explicitly; otherwise the current pipeline Python is used.
