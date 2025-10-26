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

## Database and Paths

On first execution, MimicNeoAI automatically downloads and configures the required reference database from the official FTP server.
The file `paths.yaml` is generated automatically and does not require manual editing (unless using a custom reference).

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
  pool_size: 1
  threads: 30                # base fallback
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
python -m mimicneoai.cryptic_pipeline.cryptic \
  -c /path/to/cryptic_configure.yaml
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


