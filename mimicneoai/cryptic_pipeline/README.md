# Cryptic (sORF-Encoded) Antigen Pipeline

Discover and prioritize sORF-encoded peptides from known and novel transcripts.

## Overview
1. QC
2. STAR alignment
3. Known/novel lncRNA-sORF discovery
4. Tumor/control quantification
5. HLA typing
6. Aberrantly expressed sORF peptide extraction
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
├── 07-hla_binding_pred
└── 023-shared
```

Notable subfolders:
- `04-salmon_quant/salmon_index`, `salmon_quant`, `salmon_quant_control`
- `07-hla_binding_pred/<tumor_sample>/pvacbind`

## Notes

- Tumor/control samples should be provided as `Tumor,Control` in `samples`.
- The pipeline is resumable; existing non-empty outputs are skipped.
