
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

## Database and Paths

On first execution, MimicNeoAI automatically downloads and configures the required reference database from the official FTP server.
The file `paths.yaml` is generated automatically and does not require manual editing (unless using a custom reference).

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
python -m mimicneoai.microbial_pipeline.microbial \
  -c /path/to/microbial_configure.yaml
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