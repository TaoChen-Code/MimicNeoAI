# coding=utf-8
# Script for launching neoantigen pipeline including mutation calling, annotation, HLA typing, and pVACseq processing
import os
import json
import yaml
import sys 
import argparse
sys.path.append("../functions/")
from pipline_tools import tools
from fastp import fastp
sys.path.append("./scripts/")
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool,Manager
from wgs_mouse import wgs_mouse_start
from wgs_human import wgs_human_start
from wes_mouse import wes_mouse_start
from wes_human import wes_human_start
from rna_mouse import rna_mouse_start
from annotation import annotation_vcf
from hlahd_faster import hlahd
from pvacseq import pvacseq,get_hlahd_results

# Pipeline identifier
flag = 'Neoantigen'
# Initialize utility class
log_lock = multiprocessing.Manager().Lock()  # Create a Lock object
tool = tools(os.getcwd(),flag,log_lock)
tool.write_log(f"work_path: {tool.sys_path}","info")
tool.write_log(f"start_log: {tool.start_log}","info")
tool.write_log(f"cmd_log: {tool.cmd_log}","info")

## Load configuration parameters
# Parse command line arguments
parser = argparse.ArgumentParser(description='Neo')
parser.add_argument('-c','--configure',help="Configuration file")
args = parser.parse_args()
# Load configuration
configure = tool.get_configure(args.configure)
tool.write_log(f"configures: {configure}","info")

# Load path configurations
with open(f"./configures/{flag}_pathes.yaml","r") as file:
  pathes = yaml.safe_load(file)
tool.write_log(f"pathes: {pathes}","info")

# Configure execution parameters
pool_size = configure['args']['pool_size']
species = configure['others']['species']  # human/mouse
seq_type = configure['others']['seq_type']  # wgs/wes/rna
QC = configure['others']['QC']
host_variants_calling_and_annotation = configure['others']['host_variants_calling_and_annotation']
hlatyping = configure['others']['hlatyping']
peptides_identification = configure['others']['peptides_identification']

# Process sample list
samples = configure['samples']
samples = [str(sample) for sample in samples]
tool.write_log(f"samples: {samples}","info")
tool.samples = samples  # Assign samples to utility class

# Initialize function dictionary
funcs = {"wgs_mouse":wgs_mouse_start,
         "wgs_human":wgs_human_start,
         "wes_mouse":wes_mouse_start,
         "wes_human":wes_human_start,
         "rna_mouse":rna_mouse_start}

def step1(sample,configure,pathes,tool):
    """Execute mutation calling and annotation steps"""
    # Perform mutation calling
    func = funcs[f"{seq_type}_{species}"]
    if func:
        func(sample,tool,configure,pathes)
    else:
        print("Function does not exist.")
    
    # Perform variant annotation
    if seq_type != 'rna':
        output_vep = annotation_vcf(sample,tool,pathes,configure)

def start(sample):
    """Main processing workflow for each sample"""
    sample = str(sample)
    # Quality control processing
    if QC:
        fastp(sample,configure,pathes,tool)
        
    # Variant processing pipeline
    if host_variants_calling_and_annotation:
        step1(sample,configure,pathes,tool)

    # HLA typing
    if hlatyping:
        hlahd(sample,configure,pathes,tool)    
    
    # Neoantigen prediction
    if seq_type != 'rna':
        output_dir = configure['path']['output_dir']
        output_vep = output_dir + f'/{sample}/08-vep/'
        output_hla = output_dir + f'/{sample}/09-hla/'
        pvacseq(sample,output_vep,output_hla,configure,pathes,tool)

def neoantigen():
    """Main parallel execution controller"""
    pool = NoDaemonPool(pool_size)
    for count,sample in enumerate(samples):
        pool.apply_async(start,(sample,),error_callback=tool.print_pool_error)
    pool.close()
    pool.join()
  
if __name__ == '__main__':
  tool.sharing_variable(Manager(),samples)
  neoantigen()
  tool.summary()