# coding=utf-8
# Human WES processing pipeline

"""
Pipeline wrapper for alignment and variant calling.

Required software versions (for reproducibility/documentation):
- bwa >= 0.7.17
- samtools >= 1.5
- GATK >= 4.6.0.0
- Java 17
"""

import os


def mutation_calling(run_sample_id: str, sample_name: str, tool, configure: dict, paths: dict):
    """
    Align paired-end FASTQs to the reference genome and perform preprocessing:
      1) BWA-MEM alignment -> BAM
      2) Sort BAM
      3) Mark duplicates (GATK MarkDuplicatesSpark)
      4) Base Quality Score Recalibration (BQSR)

    Args:
        run_sample_id: An identifier for logging/execution context.
        sample_name: Tumor (or Normal) sample name.
        tool: Execution helper with `judge_then_exec`, `exec_cmd`, `write_log`.
        configure: Pipeline configuration dict (paths/args/others/step_name).
        paths: Paths dict with references and tool jars.

    Returns:
        (dir_pre_vcf, dir_VQSR, bqsr_bam_path)
    """
    # --- References & tools ---
    ref = paths['database']['neoantigen']['HG38']['REF_FASTA']
    bundle_path = paths['database']['neoantigen']['HG38']['BUNDLE'] + "/"
    gatk_jar = paths['path']['common']['GATK_JAR']

    # --- I/O & runtime config ---
    tmp_dir = configure['path']['tmp_dir'] + "/"
    input_dir = configure['path']['input_dir'] + "/"
    output_dir = configure['path']['output_dir'] + "/"
    bed_file = configure['others']['bed_file']
    QC = configure['others']['QC']

    # --- Resources ---
    thread = configure['args']['thread']
    mem = configure['args']['mem']
    gatk = f"java -Xmx{mem} -jar {gatk_jar}"

    # --- Input FASTQs (QC or raw) ---
    output_dir = output_dir + sample_name + '/'
    dir_output = output_dir
    if QC:
        output_qc = dir_output + configure['step_name']['QC'] + '/'
        input_file1 = output_qc + sample_name + f'/{sample_name}.QC.R1.fq.gz'
        input_file2 = output_qc + sample_name + f'/{sample_name}.QC.R2.fq.gz'
    else:
        input_file1 = input_dir + sample_name + f'/{sample_name}.R1.fq.gz'
        input_file2 = input_dir + sample_name + f'/{sample_name}.R2.fq.gz'

    # --- Create output directory tree ---
    tool.judge_then_exec(run_sample_id, f"mkdir -p {output_dir}", output_dir)

    dir_alignment = dir_output + configure['step_name']['alignment'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_alignment}", dir_alignment)

    dir_sort = dir_output + configure['step_name']['sort'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_sort}", dir_sort)

    dir_rmdup = dir_output + configure['step_name']['rmdup'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_rmdup}", dir_rmdup)

    dir_BQSR = dir_output + configure['step_name']['bqsr'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_BQSR}", dir_BQSR)

    dir_pre_vcf = dir_output + configure['step_name']['variants_calling'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_pre_vcf}", dir_pre_vcf)

    dir_VQSR = dir_output + configure['step_name']['vqsr'] + '/'
    tool.judge_then_exec(run_sample_id, f"mkdir {dir_VQSR}", dir_VQSR)

    # --- 4.1 Alignment (BWA-MEM -> BAM) ---
    file_alignment = dir_alignment + sample_name + ".bam"
    cmd_alignment = (
        "bwa mem -t {thread} -R '@RG\\tID:foo_lane\\tPL:UNKNOWN\\tLB:library\\tSM:{sm}' {ref} {fq1} {fq2} "
        "| samtools view -S -b - > {bam}"
    ).format(thread=thread, sm=sample_name, ref=ref, fq1=input_file1, fq2=input_file2, bam=file_alignment)
    tool.judge_then_exec(run_sample_id, cmd_alignment, file_alignment)

    # --- 4.2 Sort BAM ---
    sorted_bam_filename = dir_sort + sample_name + ".sorted.bam"
    cmd_sort = f"samtools sort -@ {thread} -m 4G -O bam -o {sorted_bam_filename} {file_alignment}"
    tool.judge_then_exec(run_sample_id, cmd_sort, sorted_bam_filename)

    index_sorted_bam_filename = sorted_bam_filename + ".bai"
    cmd_sort_index = f"samtools index -@ {thread} {sorted_bam_filename}"
    tool.judge_then_exec(run_sample_id, cmd_sort_index, index_sorted_bam_filename)
    # Note: samtools -m specifies max memory per thread.

    # --- 4.3 Mark duplicates ---
    rmdup_sorted_bam_filename = dir_rmdup + sample_name + ".sorted.rmdup.bam"
    metrics_filename = dir_rmdup + sample_name + ".rmdup_metrics.txt"
    cmd_rmdup = (
        f"{gatk} MarkDuplicatesSpark --tmp-dir {tmp_dir} --spark-master local[{thread}] "
        f"-I {sorted_bam_filename} -O {rmdup_sorted_bam_filename} -M {metrics_filename}"
    )
    tool.judge_then_exec(run_sample_id, cmd_rmdup, rmdup_sorted_bam_filename)

    index_rmdup_sorted_bam_filename = rmdup_sorted_bam_filename + ".bai"
    cmd_bai = f"samtools index -@ {thread} {rmdup_sorted_bam_filename}"
    tool.judge_then_exec(run_sample_id, cmd_bai, index_rmdup_sorted_bam_filename)

    # --- 4.4 BQSR ---
    # Step 1: create recalibration table
    recal_table = dir_BQSR + sample_name + ".sorted.rmdup.recal_data.table"
    bundle_file1 = bundle_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    bundle_file2 = bundle_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    bundle_file3 = bundle_path + "dbsnp_138.hg38.vcf.gz"
    cmd_recal_table = (
        f"{gatk} BaseRecalibratorSpark --tmp-dir {tmp_dir} --spark-master local[{thread}] "
        f"-R {ref} -I {rmdup_sorted_bam_filename} "
        f"--known-sites {bundle_file1} --known-sites {bundle_file2} --known-sites {bundle_file3} "
        f"-O {recal_table} -L {bed_file}"
    )
    tool.judge_then_exec(run_sample_id, cmd_recal_table, recal_table)

    # Step 2: apply recalibration
    bqsr_bam = dir_BQSR + sample_name + ".sorted.rmdup.BQSR.bam"
    cmd_BQSR = (
        f"{gatk} ApplyBQSRSpark --tmp-dir {tmp_dir} --spark-master local[{thread}] "
        f"-R {ref} -I {rmdup_sorted_bam_filename} --bqsr-recal-file {recal_table} -O {bqsr_bam}"
    )
    tool.judge_then_exec(run_sample_id, cmd_BQSR, bqsr_bam)

    return dir_pre_vcf, dir_VQSR, bqsr_bam


def Mutect2(run_sample_id: str, sample_name: str, tool, configure: dict, paths: dict):
    """
    Run GATK Mutect2 for somatic variant calling.

    If `configure['others']['tumor_with_matched_normal']` is True,
    `sample_name` must be "TUMOR,NORMAL" (comma-separated). The function will:
      - Run alignment/preprocessing for both tumor and normal
      - Call Mutect2 in tumor-normal mode
      - Filter calls with FilterMutectCalls

    Args:
        run_sample_id: Logging/execution identifier.
        sample_name: Tumor-only name, or "tumor,normal" when matched-normal mode.
        tool: Execution helper.
        configure: Pipeline configuration.
        paths: Paths dictionary.
    """
    ref = paths['database']['neoantigen']['HG38']['REF_FASTA']
    tmp_dir = configure['path']['tmp_dir'] + "/"
    gatk_jar = paths['path']['common']['GATK_JAR']
    mem = configure['args']['mem']
    bed_file = configure['others']['bed_file']
    gatk = f"java -Xmx{mem} -jar {gatk_jar}"

    tumor_with_matched_normal = configure['others']['tumor_with_matched_normal']
    if tumor_with_matched_normal:
        tumor_sample = sample_name.split(",")[0]
        normal_sample = sample_name.split(",")[1]

        # Preprocess tumor & normal independently
        dir_pre_vcf_tumor, dir_VQSR_tumor, tumor_bam = mutation_calling(
            run_sample_id, tumor_sample, tool, configure, paths
        )
        dir_pre_vcf_normal, dir_VQSR_normal, normal_bam = mutation_calling(
            run_sample_id, normal_sample, tool, configure, paths
        )

        # Tumor + matched normal mode
        mutect_file = dir_pre_vcf_tumor + tumor_sample + '.mutect.vcf.gz'
        cmd_mutect = (
            f"{gatk} Mutect2 --tmp-dir {tmp_dir} -R {ref} "
            f"-I {tumor_bam} -I {normal_bam} -normal {normal_sample} "
            f"-O {mutect_file} -L {bed_file}"
        )
        tool.judge_then_exec(run_sample_id, cmd_mutect, mutect_file)

        filtered_mutect_file = dir_VQSR_tumor + tumor_sample + '.mutect.filtered.vcf.gz'
        cmd_mutect_filter = (
            f"{gatk} FilterMutectCalls --tmp-dir {tmp_dir} -R {ref} "
            f"-V {mutect_file} -O {filtered_mutect_file}"
        )
        tool.judge_then_exec(run_sample_id, cmd_mutect_filter, filtered_mutect_file)

    else:
        # Single-sample tumor mode is not implemented in this wrapper
        tool.write_log("Unsupported tumor_without_matched_normal", "error")
        return


def variants_calling_start(sample_name: str, tool, configure: dict, paths: dict):
    """
    Entry point to select and run a variant caller based on configuration.

    Args:
        sample_name: Sample name (or "tumor,normal" for matched mode).
        tool: Execution helper.
        configure: Pipeline configuration dict.
        paths: Paths dict.

    Notes:
        - Currently supports only Mutect2.
    """
    mutation_calling_tool = configure['others']['mutation_calling_tool']

    if mutation_calling_tool == 'Mutect2':
        Mutect2(sample_name, sample_name, tool, configure, paths)
    else:
        tool.write_log("Unsupported mutation calling tool", "error")
        return
