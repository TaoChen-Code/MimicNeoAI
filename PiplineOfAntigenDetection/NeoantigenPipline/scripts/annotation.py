# coding=utf-8
# Process VCF files from mutation calling, annotate with VEP, and identify antigen peptides
'''
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Pipeline for VCF annotation and neoantigen prediction.\n" 
    "Required software versions:\n"
    "singularity 3.9.0-rc.3\n"
    "vep110.sif\n"
    "vep102.sif\n")
'''
import os
import argparse

def annotation_vcf(sample,tool,pathes,configure):
    # Load parameters
    mouse_vep_sif = pathes['path']['mouse_vep_sif']
    human_vep_sif = pathes['path']['human_vep_sif']
    vep_data_dir = pathes['path']['vep_data_dir']
    GRCm38_ref_dir = pathes['path']['GRCm38_ref_dir']  # Mouse reference path
    hg38_ref_dir = pathes['path']['hg38_ref_dir']      # Human reference path
    species = configure['others']['species']
    output_dir = configure['path']['output_dir'] + "/"
    
    # Create output directory
    vep_dir = output_dir + '/{}/08-vep/'.format(sample)
    cmd_mkdir = "mkdir {}".format(vep_dir)
    tool.judge_then_exec(sample,cmd_mkdir,vep_dir)
    
    # Configure file paths
    mutation_calling_tool = configure['others']['mutation_calling_tool']
    if mutation_calling_tool == 'HaplotypeCaller':
        data_dir = output_dir + '/{}/07-vcf/'.format(sample)
        suffix = "HC.snps.indels.filtered.PASS" if species == "mouse" else "HC.snps.indels.VQSR.PASS"
    elif mutation_calling_tool == 'Mutect2':
        data_dir = output_dir + '/{}/06-VQSR/'.format(sample)
        suffix = "mutect.filtered"
    
    # Build absolute paths
    input_file = f"/data_dir/{sample}.{suffix}.vcf.gz"
    output_file = f"/output_dir/{sample}.{suffix}.VEP.vcf"
    abs_output_file = vep_dir + sample + f".{suffix}.VEP.vcf"
    filtered_vcf = vep_dir + sample + f".{suffix}.VEP.filtered.vcf"

    # Construct VEP command
    if species == 'mouse':
        plugins_version = 'VEP_plugins-release-102'
        ref_name = 'GRCm38.fa'
        cmd1 = f"singularity exec -B {vep_dir}:/output_dir/ -B {data_dir}:/data_dir/ -B {vep_data_dir}:/vep_data/ -B {GRCm38_ref_dir}:/fasta_data/ {mouse_vep_sif} vep --assembly GRCm38 --species mouse --input_file {input_file} --output_file {output_file} --offline --cache --dir_cache /vep_data/ --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs --fasta /fasta_data/{ref_name} --plugin Frameshift --plugin Wildtype --dir_plugins /vep_data/{plugins_version}/ --af --af_1kg --sift b"
    else:
        plugins_version = 'VEP_plugins-release-110'
        ref_name = 'hg38.fa'
        cmd1 = f"singularity exec -B {vep_dir}:/output_dir/ -B {data_dir}:/data_dir/ -B {vep_data_dir}:/vep_data/ -B {hg38_ref_dir}:/fasta_data/ {human_vep_sif} vep --input_file {input_file} --output_file {output_file} --offline --cache --dir_cache /vep_data/ --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs --fasta /fasta_data/{ref_name} --plugin Frameshift --plugin Wildtype --dir_plugins /vep_data/{plugins_version}/ --af --af_1kg --sift b"
    
    # Execute annotation
    tool.judge_then_exec(sample,cmd1,abs_output_file)
    
    # Filter mismatches
    cmd2 = f"ref-transcript-mismatch-reporter {abs_output_file} -f hard"
    tool.judge_then_exec(sample,cmd2,filtered_vcf)
    
    return vep_dir