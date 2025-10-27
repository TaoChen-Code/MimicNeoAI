
### `mimicneoai/microbial_pipeline/README.md`

```markdown
# Microbial Antigen Pipeline

Identify microbial peptides expressed in host sequencing data.

## Overview
1. Quality control (QC)  
2. Host read removal  
3. Microbial abundance quantification  
4. Microbial peptide identification  
5. HLA typing (optional)  
6. Binding and immunogenicity prediction

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
| Java (JRE/JDK)| 17              | For Java-based tools (e.g., GATK)     | `java -version`              |
| blastx       | 2.15.0+          | Microbial peptide identification      | `blastx -version`            |
| bowtie2      | 2.4.1            | Host/microbial alignment              | `bowtie2 --version`          |
| HLA-HD       | 1.7.0            | HLA typing                            | `hlahd.sh --version`         |
| apptainer    | 1.4.2            | Containerized tools (e.g., pVACtools) | `apptainer --version`        |

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

Edit only your pipeline configuration file `microbial_configure.yaml`.
An example is available at `mimicneoai/configures/microbial_configure.yaml`.

```yaml
# input_dir/
# ├── Sample_1/
# │   ├── Sample_1.R1.fq.gz
# │   └── Sample_1.R2.fq.gz
# ├── Sample_2/
# │   ├── Sample_2.R1.fq.gz
# │   └── Sample_2.R2.fq.gz
# └── Sample_3/
#     ├── Sample_3.R1.fq.gz
#     └── Sample_3.R2.fq.gz
# Configuration Paths
path:
  tmp_dir : "<project_root>/tmp/"  # Temporary directory for processing files
  input_dir : "<project_root>/input_data/"  # Input data directory
  output_dir : "<project_root>/analysis_output/"  # Output directory

# Runtime Parameters
args:
  thread : 20  # CPU threads allocated per sample
  hla_binding_threads : 5 # Number of threads for parallel pvactools runs; too many may reduce efficiency—adjust based on server performance.
  pool_size : 2  # Number of samples processed concurrently
  mem : "128G"  # Maximum memory allocation
  mem_perthread : "256M"  # Memory per thread in samtools sort

# Analysis Parameters
others:
  species : "human"  # Species selection (human)
  seq_type : "rna"  # Sequencing type (rna)
  pair : True  # Paired-end sequencing status
  QC : True  # Enable quality control steps
  host_sequences_removing : True  # Enable host sequence removal
  microbial_taxas_quantification : True  # Enable microbial taxonomy analysis
  match_length_threshold : 0.95  # Alignment length threshold
  score_threshold : 1  # Minimum alignment score
  min_pident_length : 100 # Minimum blastx match length
  microbial_peptides_identification : True  # Enable peptide identification
  hlatyping : True  # Enable HLA typing
  microbial_peptides_bindingPrediction : True  # Enable binding prediction
  # Algorithms for pVACbind (space-separated)
  algo: "BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL PickPocket SMM SMMPMBEC"

# Sample List (Note: Numeric sample IDs require quotes)
samples:
- Sample_1
- Sample_2
- Sample_3
```

## Run

```bash
# Method 1: Standard Python module execution
python -m mimicneoai.microbial_pipeline.microbial -c /path/to/microbial_configure.yaml

# Method 2: Unified CLI command (equivalent to the above)
mimicneoai microbial -c /path/to/microbial_configure.yaml
```

## Output Structure

```
<output_dir>/<sample>/
├── 00.QC
├── 01.HostSequencesRemovingStep1
├── 02.HostSequencesRemovingStep2
├── 03.MicrobialTaxasQuantificationStep1
├── 04.MicrobialTaxasQuantificationStep2
├── 05.MicrobialPeptidesIdentification
├── 06.HlaTyping
└── 07.MicrobialPeptidesBindingPrediction
```

## Notes

* Official microbial databases (taxonomy, genome, and protein indices) are downloaded and linked automatically.
* The pipeline is resumable; delete incomplete step directories before rerun if needed.