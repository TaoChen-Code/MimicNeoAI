
### `mimicneoai/mutation_derived_pipeline/README.md`

```markdown
# Mutation-Derived Neoantigen Pipeline

Detect and prioritize mutation-derived neoepitopes from WES/WGS/RNA data.

## Overview
1. Quality control (QC)  
2. Alignment and duplicate removal  
3. Base quality score recalibration (BQSR)  
4. Somatic variant calling  
5. Variant annotation (VEP)  
6. HLA typing 
7. Binding and immunogenicity prediction

## Installation
```bash
pip install mimicneoai
# or
conda install -c conda-forge -c bioconda mimicneoai
````
## External Dependencies (must be available in `PATH`)

Make sure the following tools are installed and discoverable via your `PATH` with the specified versions:

| Tool         | Required version | Purpose                               | Check command                |
|--------------|------------------|---------------------------------------|------------------------------|
| fastp        | v0.22.0          | FASTQ quality control                 | `fastp --version`            |
| bwa          | v0.7.17          | Short-read alignment                  | `bwa 2>&1 | head -n1`        |
| samtools     | v1.5             | BAM/CRAM processing                   | `samtools --version`         |
| Java (JRE/JDK)| 17              | GATK and other Java tools      | `java -version`              |
| apptainer    | 1.4.2            | Containerized tools (e.g., pVACtools) | `apptainer --version`        |
| bowtie2                 | 2.4.1            | Short-read mapping                     | `bowtie2 --version`               |
| hlahd.sh       | 1.7.0            | HLA typing                            | `hlahd.sh --version`         |

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

Edit only your pipeline configuration file `mutation_derived_configure.yaml`.
An example is available at `mimicneoai/configures/mutation_derived_configure.yaml`.

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
  thread: 30                # CPU threads per sample
  hla_binding_threads : 5 # Number of threads for parallel pvactools runs; too many may reduce efficiency—adjust based on server performance.
  pool_size: 1              # Concurrent sample processing count
  mem: "128G"               # Max memory per sample (pool_size × mem < total memory)

# Analysis Parameters
others:
  seq_type: "wes"           # Sequencing type (wes/rna)
  species: "human"          # Species (human/mouse)
  pair: True                # Paired-end sequencing (True/False)
  QC: True                  # Enable quality control

  tumor_with_matched_normal: True   # In this mode, input samples must be provided as a matched Tumor-Normal pair.
                                   # Use a comma to separate sample names, e.g. - {TumorName},{NormalName}
  host_variants_calling: True  # Host variant processing
  mutation_calling_tool: "Mutect2"  # Variant caller (Mutect2 for somatic variants)
  pad_bp: 100
  min_base_quality: 20
  min_allele_fraction: 0.02

  annotation: True

  hlatyping: True          # HLA typing enable
  peptides_identification_and_binding_prediction: True  # Peptide prediction
  mhc_i_epitope_lengths: "8,9,10,11"
  mhc_ii_epitope_lengths: "15"
  
  # WES-specific parameters
  bed_file: "/path/to/target_regions/exome_targets.bed"  # Capture regions file

  # Algorithms for pVACbind (space-separated)
  algo: "BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC"

# Sample List (Example entries)
samples:
  - Sample_001_T,Sample_001_N
  - Sample_002_T,Sample_002_N
```

## Run

```bash
# Method 1: Standard Python module execution
python -m mimicneoai.mutation_derived_pipeline.mutation_derived -c /path/to/mutation_derived_configure.yaml

# Method 2: Unified CLI command (equivalent to the above)
mimicneoai mutation-derived -c /path/to/mutation_derived_configure.yaml
```

## Output Structure

```
<output_dir>/<sample>/
├─ 00.QC/
├─ 01.alignment/
├─ 02.markdup/
├─ 03.bqsr/
├─ 04.variants_calling/
├─ 05.annotation/
├─ 06.hlatyping/
├─ 07.binding_prediction/
```
## Notes

* Tool paths, reference genomes, and database versions are managed automatically.
* Reruns skip completed steps; remove an incomplete step’s folder before rerun to ensure integrity.

