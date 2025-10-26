
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
6. HLA typing (optional)  
7. Binding and immunogenicity prediction (pVACseq)

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
  host_variants_calling_and_annotation: True  # Host variant processing
  mutation_calling_tool: "Mutect2"  # Variant caller (Mutect2 for somatic variants)
  hlatyping: False          # HLA typing enable
  peptides_identification_and_binding_prediction: False  # Peptide prediction
  
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
python -m mimicneoai.mutation_derived_pipeline.mutation-derived \
  -c /path/to/mutation_derived_configure.yaml
```

## Output Structure

```
<output_dir>/<sample>/
├─ 00.QC/
├─ 01.alignment/
├─ 02.sort/
├─ 03.rmdup/
├─ 04.BQSR/
├─ 05.variants_calling/
├─ 06.vqsr/
├─ 07.vep/
├─ 08.hlatyping/
└─ 09.binding_prediction/
```

## Notes

* Tool paths, reference genomes, and database versions are managed automatically.
* Reruns skip completed steps; remove an incomplete step’s folder before rerun to ensure integrity.

