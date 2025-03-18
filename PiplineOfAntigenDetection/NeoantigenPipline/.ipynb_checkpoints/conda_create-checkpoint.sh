conda create -n Neo 
conda activate Neo
conda install bioconda::gatk4=4.3.0.0
conda install conda-forge::openjdk=8
conda install bioconda::samtools=1.5
conda install bioconda::bwa=0.7.17
conda install bioconda::bowtie2=2.2.8
conda install conda-forge::singularity
conda install conda-forge::squashfs-tools
pip install vatools

conda create -n pvactools python=3.12
pip install pvactools
