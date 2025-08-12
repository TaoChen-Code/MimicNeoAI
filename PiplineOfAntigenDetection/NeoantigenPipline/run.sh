export TF_CPP_MIN_LOG_LEVEL=2. 
export CUDA_VISIBLE_DEVICES="".
export LD_LIBRARY_PATH=/usr/local/cuda-10.1/lib64
python Neoantigen.py -c ../../../configures/Neoantigen_configure_JLU3.yaml -p ../../../configures/Neoantigen_pathes.yaml
