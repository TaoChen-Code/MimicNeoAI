# coding=<encoding name> ： # coding=utf-8
# 处理mutation_calling得到的vcf文件，用VEP对vcf文件进行注释，然后调用pvacseq找抗原肽
#input：
#1、vcf文件s
#参数
'''
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="A simple pipeline wrapper for annotating vcf file,.\n" 
    "The software version required for this script to run: \n "
    "singularity 3.9.0-rc.3\n"
    "vep110.sif
    "vep102.sif "\n")
'''
import os
import argparse

def annotation_vcf(data_dir,input_file,sif_file,vep_data_dir,ref_dir,output_dir):
  input_file = "/data_dir/" +  vcf_file
  output_file = "/output_dir/"  +  vcf_file[:-7]  + ".VEP.vcf"
  if not os.path.exists(args.output_dir +  vcf_file[:-7]  + ".VEP.vcf"):
    #vep
    cmd1 = "singularity exec -B {}:/output_dir/ -B {}:/data_dir/ -B {}:/vep_data/ -B {}:/fasta_data/ {} vep --assembly GRCm38 --species mouse --input_file {} --output_file {} --offline --cache --dir_cache /vep_data/ --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs --fasta /fasta_data/GRCm38.fa --plugin Frameshift --plugin Wildtype --dir_plugins /vep_data/VEP_plugins-release-102/ --af --af_1kg --sift b"\
    .format(args.output_dir,vcf_dir,args.vep_data_dir,args.ref_dir,args.sif_file,input_file,output_file,)
    #print(cmd)
    exec_cmd(cmd1)  
  else:
    print(output_file+" already exists!")
    
  filtered_vcf = args.output_dir  +  vcf_file[:-7]  + ".VEP.filtered.vcf"
  if not os.path.exists(filtered_vcf):
    #mismatch
    cmd2 = "ref-transcript-mismatch-reporter {} -f hard".format(args.output_dir +  vcf_file[:-7]  + ".VEP.vcf")
    exec_cmd(cmd2)
  else:
    print(output_file+" already exists!")
  



  




