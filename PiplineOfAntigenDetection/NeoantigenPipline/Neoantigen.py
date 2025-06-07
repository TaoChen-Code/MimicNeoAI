# coding=utf-8
# Script for launching neoantigen pipeline including mutation calling, annotation, HLA typing, and pVACseq processing
import os
import json
import yaml
import sys 
import argparse
import traceback
sys.path.append("../functions/")
from pipline_tools import tools
from fastp import fastp
sys.path.append("./scripts/")
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool,Manager
# from wgs_mouse import wgs_mouse_start
# from wgs_human import wgs_human_start
# from wes_mouse import wes_mouse_start
from wes_human import wes_human_start
# from rna_mouse import rna_mouse_start
from annotation import annotation_vcf
from hlahd_faster import hlahd
from pvacseq import Pvacseq

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)
    
# ##python3.11
# class NoDaemonPool(multiprocessing.pool.Pool):
#     @staticmethod
#     def Process(_, *args, **kwds):
#         return NoDaemonProcess(*args, **kwds)

#python3.7
class NoDaemonPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

# Pipeline identifier
flag = 'Neoantigen'
# Initialize utility class
manager = multiprocessing.Manager()
shared_lock = manager.Lock()  # This lock CAN be pickled
tool = tools(os.getcwd(),flag,shared_lock)
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
funcs = {#"wgs_mouse":wgs_mouse_start,
         #"wgs_human":wgs_human_start,
         #"wes_mouse":wes_mouse_start,
         "wes_human":wes_human_start,
         #"rna_mouse":rna_mouse_start
        }

# Configure pipeline step names
configure['step_name'] = {}
configure['step_name']['QC'] = '00.QC'
configure['step_name']['alignment'] = '01.alignment'
configure['step_name']['sort'] = '02.sort'
configure['step_name']['rmdup'] = '03.rmdup'
configure['step_name']['bqsr'] = '04.BQSR'
configure['step_name']['variants_calling'] = '05.variants_calling'
configure['step_name']['vqsr'] = '06.vqsr'
configure['step_name']['annotation'] = '07.vep'
configure['step_name']['hla'] = '08.hlatyping'
configure['step_name']['pvacseq'] = '09.peptides_identification' 

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
    try:
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
        if seq_type != 'rna' and peptides_identification:
            output_dir = configure['path']['output_dir']
            output_vep = output_dir + f'/{sample}/08-vep/'
            step_name_hla = configure['step_name']['hla']
            output_hla = output_dir + f'/{sample}/{step_name_hla}/'
            pvacseq = Pvacseq(tool)
            pvacseq.run_pvacseq_parallel(sample,output_vep,output_hla,configure,pathes)
    except Exception as e:
        error_message = traceback.format_exc()
        #print(f"start Error occurred: {error_message}")
        tool.write_log(f"start Error occurred: {error_message}","error")

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