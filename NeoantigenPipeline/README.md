# Neoantigen
## install 

### **1. Download release**, pkgs , refs and test data

#### (1) release:

####  [MimicNeoAI.zip]()

```shell
git clone https://github.com/BioStaCs-public/MimicNeoAI.git
```

#### (2) pkgs:

####   [gatk-4.6.0.0.zip](https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip)

```shell
#download pkgs
YOUR_PATH="/your_custom/"
mkdir -p "${YOUR_PATH}/pkgs"
cd "${YOUR_PATH}/pkgs"
wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
unzip gatk-4.6.0.0.zip
```
Download and install hlahd1.7.0 from https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/

Installing IEDB binding prediction tools
https://pvactools.readthedocs.io/en/latest/install.html

#### (3) refs and test data: 

https://doi.org/10.5281/zenodo.15582924

### 2. create conda env

Successed in 20.04.1-Ubuntu x86_64 GNU/Linux
```shell
conda create -n java17
conda activate java17
conda install conda-forge::openjdk=17

conda create -n Neoantigen python=3.6
pip install pyyaml tqdm pandas vatools
conda install -y bioconda::fastp=0.22.0
conda install -y bioconda::samtools=1.5
("ln -sf ~/.conda/envs/Neoantigen/lib/libcrypto.so.1.1 ~/.conda/envs/Neoantigen/lib/libcrypto.so.1.0.0" if error while loading shared libraries: libcrypto.so.1.0.0)
conda install -y bioconda::bwa=0.7.17
conda install -y bioconda::bowtie2=2.3.5.1
conda install -y conda-forge::singularity
conda install -y bioconda::bcftools
("ln -s /lib/x86_64-linux-gnu/libgsl.so.23 ~/.conda/envs/Neoantigen/lib/libgsl.so.25")
pip install pvactools==4.2.1
pip install mhcflurry
mhcflurry-downloads fetch( Download the data manually: mhcflurry-downloads fetch models_class1_presentation --already-downloaded-dir downloads)
pip install tensorflow==1.15.0
conda install -c conda-forge cudatoolkit=10.0 cudnn=7.6
pip install mhcnuggets
pip install git+https://github.com/griffithlab/bigmhc.git#egg=bigmhc
BIGMHC_DIR=$(dirname $(which bigmhc_predict))
echo "export PATH=\$PATH:$BIGMHC_DIR" >> ~/.bashrc
source ~/.bashrc
pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno
DEEPIMMUNO_DIR=$(dirname $(which deepimmuno-cnn))
echo "export PATH=\$PATH:$DEEPIMMUNO_DIR" >> ~/.bashrc
source ~/.bashrc
```

### 3. ‌**Edit the configuration file**

#### (1)  ‌‌**Use the `conda env list` command**‌ to display the file paths of our Conda environments on your system.

Execute the following command to view all created virtual environments and their paths:

```shell
conda env list
# or
conda info --envs
```

Example output:

```shell
# conda environments:
base                  /home/user/anaconda3
myenv                 /home/user/anaconda3/envs/myenv
*MicrobialAntigen   ~/.conda/envs/MicrobialAntigen/
java17                ~/.conda/envs/java17/
```

The `*` symbol in the list marks the currently activated environment, and the path column shows the locations of all virtual environments‌

#### (2) **Configure in the following files**

**Configure1：configures/Neoantigen_pathes.yaml**

(Use absolute paths if you're not sure how Python handles relative paths.)

```yaml
# Neoantigen_pathes.yaml
# Resource Path Configuration
path:
  # Human reference genome resources
  hg38_ref_dir: "/path/to/neoantigen/WGS/hg38/"
  hg38_ref_file: "/path/to/neoantigen/WGS/hg38/hg38.fa"
  hg38_bundle_path: "/path/to/neoantigen/WGS/human_bundle/"

  # Tool executables
  java17: '/opt/conda/envs/java17/bin/java'
  gatk_4.6_jar: '/opt/biotools/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar'
  fastp: '/opt/conda/envs/scMic/bin/fastp'

  # VEP resources
  vep_data_dir: "/path/to/neoantigen/pvactools/vep_data"
  human_vep_sif: "/path/to/neoantigen/pvactools/vep.sif"

  # HLA typing resources
  freq_data_dir: '/path/to/hlahd1.7.0_new/hlahd.1.7.0/freq_data/'
  HLA_gene: '/path/to/hlahd1.7.0_new/hlahd.1.7.0/HLA_gene.split.txt'
  dictionary: '/path/to/hlahd1.7.0_new/hlahd.1.7.0/dictionary/'
  picard_path: "/path/to/BioDatahub/picard.jar"
  hla_gen: "/path/to/hlahd1.7.0_new/hla_gen/hla_gen"

  # Antigen prediction configuration
  pvacseq: '/opt/conda/envs/pvactools/bin/pvacseq'
  iedb-install-directory: '/opt/db/IEDB/'

  # R
  R_HOME: '/opt/conda/envs/r4.3/lib/R/'
  R_LIBRARY: "/opt/conda/envs/r4.3/lib/R/library"
```

**Configure2：configures/Neoantigen_configure.yaml**

(Single-end data must meet the following requirements: the folder should be named after the sample (sample), and it should contain `sample.fq.gz`, with the suffix required to be `.fq.gz`. The file format should be `sample/sample.fq.gz`. For paired-end data, the file formats should be `sample/sample.R1.fq.gz` and `sample/sample.R2.fq.gz`, with the suffixes required to be `.R1.fq.gz` and `.R2.fq.gz`, respectively.)

```yaml
# Neoantigen_configure.yaml
# Configuration Paths
path:
  tmp_dir: "/path/to/project/tmp/"          # Directory for temporary files
  input_dir: "/path/to/raw_data/"           # Input data directory containing samples
  output_dir: "/path/to/analysis_results/" # Output directory for final results

# Runtime Parameters
args:
  thread: 30                # CPU threads per sample
  pool_size: 2              # Concurrent sample processing count
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
  peptides_identification: False  # Peptide prediction
  
  # WES-specific parameters
  bed_file: "/path/to/target_regions/exome_targets.bed"  # Capture regions file
  time_out: 36000          # Maximum runtime in seconds

# Sample List (Example entries)
samples:
  - Sample_001
  - Sample_002
```

### **4. Run pipline on test data**

```shell
conda activate Neoantigen
YOUR_PATH="/your_custom/"
cd "${YOUR_PATH}/MimicNeoAI/NeoantigenPipline/"
./run.sh
# or
export TF_CPP_MIN_LOG_LEVEL=2. 
export CUDA_VISIBLE_DEVICES="".
python Neoantigen.py -c ../configures/Neoantigen_configure.yaml -p ../configures/Neoantigen_pathes.yaml
```

**Runtime logs can be viewed in the command-line terminal and are backed up in `${YOUR_PATH}/MimicNeoAI/NeoantigenPipline/log`


### **5. Important Notes**‌
**① Pipeline Interruption Handling**‌
If the pipeline is terminated due to errors or manually interrupted, it will ‌skip existing non-empty files‌ upon restarting. To ensure pipeline integrity:

‌**Mandatory verification‌:**
Before rerunning, confirm all pre-existing files (e.g., in <project_root>/sample/01.alignment/) are fully generated and valid.
**‌Uncertainty protocol‌:**
If file completeness cannot be verified, delete them using:
```shell
rm -rf <project_root>/tmp/sample_01/alignment/*
```
‌**Critical warning‌:**
Residual corrupted files may cause downstream analysis failures or data contamination.
