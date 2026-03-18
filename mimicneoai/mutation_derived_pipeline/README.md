# Mutation-Derived Neoantigen Pipeline

Detect and prioritize mutation-derived neoepitopes from matched tumor-normal sequencing data.

## Overview
1. QC
2. Alignment / duplicate marking / BQSR
3. Parallel somatic calling (Mutect2 + Strelka2 + VarDict)
4. VEP annotation and filtering
5. HLA typing
6. pVACseq-based binding prediction

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
| `bwa` | `0.7.17` | Alignment |
| `samtools` | `1.5` | BAM/SAM processing |
| `java` | `17` | Run GATK JAR |
| `GATK` (jar) | `4.6.0.0` | Mutect2/MarkDuplicates/BQSR |
| `bcftools` | `1.11` | VCF filtering/normalization/indexing |
| `bgzip`, `tabix` | `htslib 1.11+` | BED/VCF compression and indexing |
| `configManta.py` | from your Manta install | Manta workflow configuration |
| `configureStrelkaSomaticWorkflow.py` | from your Strelka2 install | Strelka2 workflow configuration |
| `vardict-java` | compatible VarDict release | VarDict caller |
| `testsomatic.R` | from your VarDict install | VarDict post-processing |
| `var2vcf_paired.pl` | from your VarDict install | VarDict VCF conversion |
| `ref-transcript-mismatch-reporter` | compatible pVACtools release | Post-VEP mismatch filtering |
| `VEP` (container/plugins) | VEP image in `paths.yaml`; plugins `VEP_plugins-release-110` | Variant annotation |
| `bowtie2` | `2.4.1` | HLA typing pre-alignment |
| `hlahd.sh` | `1.7.0` | HLA typing |
| `apptainer` | `1.4.2` | VEP and pVACtools container execution |
| `pVACtools` (container) | `4.2.1` | Binding prediction (`pvacseq`) |
| `awk` | system tool | FASTQ header normalization in HLA typing |

## Database and Paths

Shared documentation:
- [`mimicneoai/configures/Database_and_Paths.md`](../configures/Database_and_Paths.md)

## Configuration

Use the canonical template directly:
- [`mimicneoai/configures/mutation_derived_configure.yaml`](../configures/mutation_derived_configure.yaml)

Important behavior from current code:
- `others.tumor_with_matched_normal` must be `True`.
- `samples` should be `Tumor,Normal` pairs.

## Run

```bash
# Method 1
python -m mimicneoai.mutation_derived_pipeline.mutation_derived -c /path/to/mutation_derived_configure.yaml

# Method 2 (unified CLI)
mimicneoai mutation-derived -c /path/to/mutation_derived_configure.yaml
```

## Output Structure

Pipeline outputs are written under:

```text
<output_dir>/Mutation-derived/<tumor_sample>/
├── 00.QC
├── 01.alignment
├── 02.markdup
├── 03.bqsr
├── 04.variants_calling
├── 05.annotation
├── 06.hlatyping
└── 07.binding_prediction
```

Notable subfolders:
- `04.variants_calling/Mutect2`
- `04.variants_calling/Strelka2`
- `04.variants_calling/VarDict`
- `04.variants_calling/Merge3Callers`
- `05.annotation/chunks`
- `07.binding_prediction/chunks` and `07.binding_prediction/combined`

## Notes

- This pipeline currently runs matched tumor-normal mode only.
- The pipeline is resumable; existing non-empty outputs are skipped.
