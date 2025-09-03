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
import gzip
from typing import Union
import argparse


def split_vcf(input_vcf: str, output_dir: str, file_suffix: str, chunks: int = 10) -> None:
    """
    Split a VCF file into multiple chunks with equal distribution of variants.

    Parameters:
        input_vcf (str): Path to the input VCF file (can be .vcf or .vcf.gz)
        chunks (int): Number of chunks to split into (default: 10)

    Returns:
        None: Writes output files in the current directory

    Output Files:
        chunk_{0-N}.vcf or chunk_{0-N}.vcf.gz (matches input format)
    """

    # Step 1: Read header lines and variant lines separately
    header = []
    variants = []

    # Open file with appropriate handler (gzip or regular)
    file_handler = gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf)

    with file_handler as f:
        for line in f:
            if line.startswith('#'):  # Header line
                header.append(line)
            else:  # Variant line
                variants.append(line)

    # Step 2: Calculate distribution of variants across chunks
    total_variants = len(variants)
    base_lines_per_chunk = total_variants // chunks  # Minimum variants per chunk
    remainder = total_variants % chunks  # Leftover variants to distribute

    # Step 3: Create chunks with distributed variants
    chunk_start = 0
    for i in range(chunks):
        # Determine this chunk's size (base + 1 if we have remainders left)
        current_chunk_size = base_lines_per_chunk + (1 if i < remainder else 0)
        chunk_end = chunk_start + current_chunk_size

        # Extract variants for this chunk
        chunk_variants = variants[chunk_start:chunk_end]

        # Determine output format (matches input)
        output_ext = '.vcf.gz' if input_vcf.endswith('.gz') else '.vcf'
        output_file = f"{output_dir}/{file_suffix}.chunk_{i}{output_ext}"

        # Write output with appropriate handler
        writer = gzip.open(output_file, 'wt') if output_ext == '.vcf.gz' else open(output_file, 'w')
        with writer as f:
            f.writelines(header)  # Write header lines
            f.writelines(chunk_variants)  # Write variant lines

        # Update pointer for next chunk
        chunk_start = chunk_end


def annotation_vcf(run_sample_id, sample, tool, pathes, configure):
    # Load parameters
    mouse_vep_sif = pathes['path']['mouse_vep_sif']
    human_vep_sif = pathes['path']['human_vep_sif']
    vep_data_dir = pathes['path']['vep_data_dir']
    GRCm38_ref_dir = pathes['path']['GRCm38_ref_dir']  # Mouse reference path
    hg38_ref_dir = pathes['path']['hg38_ref_dir']  # Human reference path
    species = configure['others']['species']
    output_dir = configure['path']['output_dir'] + "/"
    thread = configure['args']['thread']

    # Create output directory
    vep_dir = output_dir + '/{}/'.format(sample) + configure['step_name']['annotation'] + '/'
    cmd_mkdir = "mkdir {}".format(vep_dir)
    tool.judge_then_exec(run_sample_id, cmd_mkdir, vep_dir)

    # Configure file paths
    mutation_calling_tool = configure['others']['mutation_calling_tool']
    data_dir = output_dir + '/{}/'.format(sample) + configure['step_name']['vqsr'] + '/'
    if mutation_calling_tool == 'HaplotypeCaller':
        suffix = "HC.snps.indels.filtered.PASS" if species == "mouse" else "HC.snps.indels.VQSR.PASS"
    elif mutation_calling_tool == 'Mutect2':
        suffix = "mutect.filtered"

    # Build absolute paths
    input_file = f"/data_dir/{sample}.{suffix}.vcf.gz"
    output_file = f"/output_dir/{sample}.{suffix}.VEP.vcf"
    abs_output_file = vep_dir + sample + f".{suffix}.VEP.vcf"
    filtered_vcf = f"/output_dir/{sample}.{suffix}.VEP.filtered.vcf"
    abs_filtered_vcf = vep_dir + sample + f".{suffix}.VEP.filtered.vcf"
    output_file_tsv = f"/output_dir/{sample}.{suffix}.VEP.filtered.tsv"
    abs_output_file_tsv = vep_dir + sample + f".{suffix}.VEP.filtered.tsv"

    '''
    The --pick option might be useful to limit the annotation to the “top” transcript for each variant (the one for which the most dramatic consequence is predicted). Otherwise, VEP will annotate each variant with all possible transcripts. pVACseq will provide predictions for all transcripts in the VEP CSQ field. Running VEP without the --pick option can therefore drastically increase the runtime of pVACseq.
    The --transcript_version option will add the transcript version to the transcript identifiers. This option might be needed if you intend to annotate your VCF with expression information. Particularly if your expression estimation tool uses versioned transcript identifiers (e.g. ENST00000256474.2).
    '''
    # Construct VEP command
    vep_options = f"--offline --cache --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs --plugin Frameshift --plugin Wildtype --af --af_1kg --sift b --transcript_version --pick"
    if species == 'mouse':
        plugins_version = 'VEP_plugins-release-102'
        ref_name = 'GRCm38.fa'
        cmd1 = f"singularity exec -B {vep_dir}:/output_dir/ -B {data_dir}:/data_dir/ -B {vep_data_dir}:/vep_data/ -B {GRCm38_ref_dir}:/fasta_data/ {mouse_vep_sif} vep --fork {thread} --assembly GRCm38 --species mouse --input_file {input_file} --output_file {output_file} --dir_cache /vep_data/ --fasta /fasta_data/{ref_name}  --dir_plugins /vep_data/{plugins_version}/ {vep_options}"
    else:
        plugins_version = 'VEP_plugins-release-110'
        ref_name = 'hg38.fa'
        cmd1 = f"singularity exec -B {vep_dir}:/output_dir/ -B {data_dir}:/data_dir/ -B {vep_data_dir}:/vep_data/ -B {hg38_ref_dir}:/fasta_data/ {human_vep_sif} vep --fork {thread} --input_file {input_file} --output_file {output_file} --dir_cache /vep_data/ --fasta /fasta_data/{ref_name} --dir_plugins /vep_data/{plugins_version}/ {vep_options}"

    # Execute annotation
    tool.judge_then_exec(run_sample_id, cmd1, abs_output_file)

    # Filter mismatches
    cmd2 = f"ref-transcript-mismatch-reporter {abs_output_file} -f hard"
    tool.judge_then_exec(run_sample_id, cmd2, abs_filtered_vcf)

    # === NEW STEP: keep only FILTER==PASS before chunking ===
    abs_filtered_pass_vcf = vep_dir + sample + f".{suffix}.VEP.filtered.PASS.vcf"
    cmd_keep_pass = (
        f"awk 'BEGIN{{FS=OFS=\"\\t\"}} /^#/ || $7==\"PASS\"' "
        f"{abs_filtered_vcf} > {abs_filtered_pass_vcf}"
    )
    tool.write_log("Filter VCF to FILTER==PASS for downstream analysis.", "info")
    tool.judge_then_exec(run_sample_id, cmd_keep_pass, abs_filtered_pass_vcf)

    # split vcf
    chunk_dir = vep_dir + '/chunks'
    cmd_mkdir = "mkdir {}".format(chunk_dir)
    tool.judge_then_exec(run_sample_id, cmd_mkdir, chunk_dir)

    tool.write_log(f"split PASS-only vcf to {thread} chunks.", "info")
    split_vcf(abs_filtered_pass_vcf, chunk_dir, f"{sample}.{suffix}.VEP.filtered.PASS", chunks=thread)
    tool.write_log("split vcf done!", "info")
    return vep_dir