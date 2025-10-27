### `mimicneoai/cryptic_pipeline/README.md`

````markdown
# Cryptic (sORF-Encoded) Antigen Pipeline

Discover and prioritize short ORF (sORF)-encoded peptides from known and novel transcripts.

## Overview
1. Quality control (QC)  
2. Alignment and transcript assembly  
3. sORF discovery  
4. Tumor/control expression quantification  
5. Aberrant (tumor-specific) expression filtering  
6. HLA typing (optional)  
7. Binding and immunogenicity prediction

## Installation
```bash
pip install mimicneoai
# or
conda install -c conda-forge -c bioconda mimicneoai
````

## External Dependencies (must be available in `PATH`)

Make sure the following tools are installed and discoverable via your `PATH` with the specified versions:

| Tool                    | Required version | Purpose                                | Check command                     |
|--------------------------|------------------|----------------------------------------|-----------------------------------|
| fastp                   | v0.22.0          | FASTQ quality control                  | `fastp --version`                 |
| STAR                    | 2.5.3a           | RNA alignment                          | `STAR --version`                  |
| samtools                | v1.5             | BAM/CRAM processing                    | `samtools --version`              |
| stringtie               | 3.0.1            | Transcript assembly                    | `stringtie --version`             |
| gffcompare              | v0.12.10         | Transcript annotation comparison       | `gffcompare -v`                   |
| gffread                 | 0.12.7           | GFF/GTF file processing                | `gffread --version`               |
| minimap2                | 2.30-r1287       | Transcript-to-reference alignment      | `minimap2 --version`              |
| bcftools                | 1.11 (htslib 1.13)| Variant calling and processing         | `bcftools --version`              |
| TransDecoder.LongOrfs   | 5.5.0            | ORF prediction                         | `TransDecoder.LongOrfs -h`        |
| salmon                  | 1.10.0           | Transcript quantification              | `salmon --version`                |
| bowtie2                 | 2.4.1            | Short-read mapping                     | `bowtie2 --version`               |
| hlahd.sh                  | 1.7.0            | HLA typing                             | `hlahd.sh --version`              |
| apptainer               | 1.4.2            | Containerized tools (e.g., pVACtools)  | `apptainer --version`             |

> Tip: If a command above is not found or the version is lower than required, install/upgrade it and ensure the binary is on your `PATH`.


## Database and Paths

On first execution, MimicNeoAI automatically downloads and configures the required reference database from the official FTP server.
 The file `paths.yaml` is generated automatically and does not require manual editing (unless using a custom reference).

### 🧩 Database Manual Download

If the automatic download fails or you prefer to manage storage manually,
 you can manually download and extract the MimicNeoAI database using the built-in helper:

#### 1. Default download to the project path (`mimicneoai/database`)

```bash
# Method 1: Standard Python module execution
python -m mimicneoai.download_database

# Method 2: Unified CLI command (equivalent to the above)
mimicneoai download_database
```

#### 2. Custom download path (recommended for limited system disk space)

```bash
# Method 1: Specify a custom directory for download and extraction
python -m mimicneoai.download_database --target-dir /mnt/data/MimicNeoAI_DB

# Method 2: Equivalent CLI command
mimicneoai download_database --target-dir /mnt/data/MimicNeoAI_DB
```

If a custom path is specified, the extracted folder will be **symlinked to `mimicneoai/database/`** automatically.

## Configuration

Edit only the pipeline configuration file `cryptic_configure.yaml`.
An example is available at `mimicneoai/configures/cryptic_configure.yaml`.

```yaml
# input_dir/
# ├── Sample_001/
# │   ├── Sample_001.R1.fq.gz
# │   └── Sample_001.R2.fq.gz
# ├── Sample_002/
# │   ├── Sample_002.R1.fq.gz
# │   └── Sample_002.R2.fq.gz
# └── Sample_003/
#     ├── Sample_003.R1.fq.gz
#     └── Sample_003.R2.fq.gz
# Configuration Paths
path:
  tmp_dir: "/path/to/project/tmp/"          # Directory for temporary files
  input_dir: "/path/to/raw_data/"           # Input data directory containing samples
  output_dir: "/path/to/analysis_results/" # Output directory for final results

# Runtime Parameters
args:
  threads: 30   # base fallback
  hla_binding_threads : 5 # Number of threads for parallel pvactools runs; too many may reduce efficiency—adjust based on server performance.
  pool_size: 1
  trinity_mem: "100G"
others:
  QC: true
  alignment: true
  known: true
  novel: true
  salmon_quant: true
  salmon_quant_control: true
  hlatyping: true
  extract_aeseps: true
  hla_binding_pred: true

  min_tpm_tumor: 5.0
  max_tpm_ctrl: 0.5
  min_log2fc: 4.0
  # Algorithms for pVACbind (space-separated)
  algo: "BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC"
  mhcI_lengths: "8,9,10"
  mhcII_lengths: "15"

samples:
  - Sample_001_T,Sample_001_N
  - Sample_002_T,Sample_002_N
```

## Run

```bash
# Method 1: Standard Python module execution
python -m mimicneoai.cryptic_pipeline.cryptic -c /path/to/cryptic_configure.yaml

# Method 2: Unified CLI command (equivalent to the above)
mimicneoai cryptic -c /path/to/cryptic_configure.yaml
```

## Output Structure

```
<output_dir>/<sample>/
├─ 00.QC/
├─ 01.alignment/
├─ 02.sORF_discovery/
├─ 03.salmon_tumor/
├─ 04.salmon_control/
├─ 05.hla_typing/
├─ 06.sORF_encoded_peptides/
└─ 07.binding_prediction/
```

## Notes

* The pipeline is resumable: existing non-empty outputs are skipped.
* If a step fails or is incomplete, delete that step’s folder and rerun.
* Tumor/control pairs must be listed as `Tumor,Control`.


