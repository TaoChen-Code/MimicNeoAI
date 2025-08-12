# coding=utf-8
# Human WES processing pipeline
'''
"Pipeline wrapper for alignment and variant calling.\n"
    "Required software versions:\n"
    "bwa-0.7.17-r1188\n"
    "samtools-1.5\n"
    "picard-2.25.5\n"
    "gatk-4.6.0.0\n"
    "java17"
'''
import os


def mutation_calling(run_sample_id,sample_name,tool,configure,pathes):
    #1 Read input parameters
    ref = pathes['path']['hg38_ref_file']
    bundle_path = pathes['path']['hg38_bundle_path']
    java17 = pathes['path']['java17']
    gatk_jar = pathes['path']['gatk_4.6_jar']
    
    # Path configurations
    tmp_dir = configure['path']['tmp_dir'] + "/"  # Temporary files directory
    input_dir = configure['path']['input_dir'] + "/"  # Input files directory
    output_dir = configure['path']['output_dir'] + "/"  # Output directory
    bed_file = configure['others']['bed_file']  # Target regions BED file
    mutation_calling_tool = configure['others']['mutation_calling_tool']
    QC = configure['others']['QC']
    
    # Resource allocation
    thread = configure['args']['thread']  # CPU threads per sample
    mem = configure['args']['mem']  # Max memory for GATK
    gatk = "{} -Xmx{} -jar {}".format(java17,mem,gatk_jar)

    #2 Initialize output structure
    output_dir = output_dir + sample_name +'/'
    dir_output = output_dir
    if QC:
        output_qc = output_dir + '/' + configure['step_name']['QC'] + '/'
        input_file1 = output_qc + sample_name + '/' + "{}.QC.R1.fq.gz".format(sample_name)
        input_file2 = output_qc + sample_name + '/' + "{}.QC.R2.fq.gz".format(sample_name)
    else:
        input_file1 = output_dir + sample_name + '/' + "{}.R1.fq.gz".format(sample_name)
        input_file2 = output_dir + sample_name + '/' + "{}.R2.fq.gz".format(sample_name)
    
    #3 Create directory structure
    ## Main output directory
    tool.judge_then_exec(run_sample_id,"mkdir -p {}".format(output_dir),output_dir)
    
    ### Create processing directories
    dir_alignment = dir_output + configure['step_name']['alignment'] + '/'
    cmd_mkdir_alignment = "mkdir {}".format(dir_alignment) 
    tool.judge_then_exec(run_sample_id,cmd_mkdir_alignment,dir_alignment)

    dir_sort = dir_output + configure['step_name']['sort'] + '/'
    cmd_mkdir_sort = "mkdir {}".format(dir_sort)
    tool.judge_then_exec(run_sample_id,cmd_mkdir_sort,dir_sort)

    dir_rmdup = dir_output + configure['step_name']['rmdup'] + '/'
    cmd_mkdir_rmdup = "mkdir {}".format(dir_rmdup)
    tool.judge_then_exec(run_sample_id,cmd_mkdir_rmdup,dir_rmdup)

    dir_BQSR = dir_output + configure['step_name']['bqsr'] +'/'
    cmd_mkdir_BQSR = "mkdir {}".format(dir_BQSR)
    tool.judge_then_exec(run_sample_id,cmd_mkdir_BQSR,dir_BQSR)

    dir_pre_vcf = dir_output + configure['step_name']['variants_calling'] + '/'
    cmd_mkdir_pre_vcf = "mkdir {}".format(dir_pre_vcf)
    tool.judge_then_exec(run_sample_id,cmd_mkdir_pre_vcf,dir_pre_vcf)

    dir_VQSR = dir_output + configure['step_name']['vqsr'] + '/'
    cmd_mkdir_VQSR = "mkdir {}".format(dir_VQSR)
    tool.judge_then_exec(run_sample_id,cmd_mkdir_VQSR,dir_VQSR)

    #4 Data processing pipeline
    ##4.1 Alignment with BWA
    file_alignment = dir_alignment+sample_name+".bam"
    cmd_alignment = "bwa mem -t {} -R '@RG\\tID:foo_lane\\tPL:UNKNOWN\\tLB:library\\tSM:{}' {} {} {} | samtools view -S -b - > {}"\
    .format(thread,sample_name,ref,input_file1,input_file2,file_alignment)
    tool.judge_then_exec(run_sample_id,cmd_alignment,file_alignment)

    ##4.2 Sort BAM file
    sorted_bam_filename = dir_sort+sample_name+".sorted.bam"
    cmd_sort = "samtools sort -@ {} -m 4G -O bam -o {} {}"\
    .format(thread,sorted_bam_filename,file_alignment)
    tool.judge_then_exec(run_sample_id,cmd_sort,sorted_bam_filename)

    index_sorted_bam_filename = dir_sort+sample_name+".sorted.bam.bai"
    cmd_sort_index = "samtools index -@ {} {}".format(thread,sorted_bam_filename)
    tool.judge_then_exec(run_sample_id,cmd_sort_index,index_sorted_bam_filename)
    # samtools -m specifies max memory per thread

    ##4.3 Remove duplicates
    rmdup_sorted_bam_filename = dir_rmdup+sample_name+".sorted.rmdup.bam"
    metrics_filename = dir_rmdup+sample_name+".rmdup_metrics.txt"
    cmd_rmdup = "{} MarkDuplicatesSpark --tmp-dir {} --spark-master local[{}] -I {} -O {} -M {}"\
    .format(gatk,tmp_dir,thread,sorted_bam_filename,rmdup_sorted_bam_filename,metrics_filename)
    tool.judge_then_exec(run_sample_id,cmd_rmdup,rmdup_sorted_bam_filename)

    index_rmdup_sorted_bam_filename = dir_rmdup+sample_name+".sorted.rmdup.bam.bai"
    cmd_bai = "samtools index -@ {} {}".format(thread,rmdup_sorted_bam_filename)
    tool.judge_then_exec(run_sample_id,cmd_bai,index_rmdup_sorted_bam_filename)

    ##4.4 Base Quality Score Recalibration (BQSR)
    ### Step1: Generate recalibration table
    recal_table_index_rmdup_sorted_bam_filename = dir_BQSR+sample_name+".sorted.rmdup.recal_data.table"
    bundle_file1 = bundle_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    bundle_file2 = bundle_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    bundle_file3 = bundle_path + "dbsnp_138.hg38.vcf.gz"
    cmd_recal_table = "{} BaseRecalibratorSpark --tmp-dir {} --spark-master local[{}] -R {} -I {} --known-sites {} --known-sites {} --known-sites {} -O {} -L {}"\
    .format(gatk,tmp_dir,thread,ref,rmdup_sorted_bam_filename,bundle_file1,bundle_file2,bundle_file3,recal_table_index_rmdup_sorted_bam_filename,bed_file)
    tool.judge_then_exec(run_sample_id,cmd_recal_table,recal_table_index_rmdup_sorted_bam_filename)

    ### Step2: Apply BQSR
    BQSR_index_rmdup_sorted_bam_filename = dir_BQSR+sample_name+".sorted.rmdup.BQSR.bam"
    cmd_BQSR = "{} ApplyBQSRSpark --tmp-dir {} --spark-master local[{}] -R {} -I {} --bqsr-recal-file {} -O {}"\
    .format(gatk,tmp_dir,thread,ref,rmdup_sorted_bam_filename,recal_table_index_rmdup_sorted_bam_filename,BQSR_index_rmdup_sorted_bam_filename)
    tool.judge_then_exec(run_sample_id,cmd_BQSR,BQSR_index_rmdup_sorted_bam_filename)

    #5 Variant calling
    ##5.1 Initial variant calling
    if mutation_calling_tool == 'HaplotypeCaller':
        ### HaplotypeCaller workflow
        gvcf_file = dir_pre_vcf + sample_name + ".g.vcf.gz"
        cmd_gvcf = "{} HaplotypeCaller --tmp-dir {} -R {} -I {} --emit-ref-confidence GVCF -O {} -L {}"\
        .format(gatk,tmp_dir,ref,BQSR_index_rmdup_sorted_bam_filename,gvcf_file,bed_file)
        tool.judge_then_exec(run_sample_id,cmd_gvcf,gvcf_file)

        vcf_file = dir_pre_vcf + sample_name + ".HC.vcf.gz"
        cmd_vcf = "{} GenotypeGVCFs --tmp-dir {} -R {} --variant {} -O {}"\
        .format(gatk,tmp_dir,ref,gvcf_file,vcf_file)
        tool.judge_then_exec(run_sample_id,cmd_vcf,vcf_file)

        ##5.2 Variant Quality Score Recalibration (VQSR)
        ### SNP Recalibration
        snps_recal_file = dir_VQSR + sample_name + ".HC.snps.recal"
        tranches_file = dir_VQSR + sample_name + ".HC.snps.tranches"
        bundle_file1 = bundle_path + "hapmap_3.3.hg38.vcf.gz"
        bundle_file2 = bundle_path + "1000G_omni2.5.hg38.vcf.gz"
        bundle_file3 = bundle_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        bundle_file4 = bundle_path + "dbsnp_138.hg38.vcf.gz"
        cmd_recal = "{} VariantRecalibrator --tmp-dir {} -R {} --variant {} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {} -resource:omini,known=false,training=true,truth=false,prior=12.0 {} -resource:1000G,known=false,training=true,truth=false,prior=10.0 {} -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 {} -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP --output {} -tranches-file {}"\
        .format(gatk,tmp_dir,ref,vcf_file,bundle_file1,bundle_file2,bundle_file3,bundle_file4,snps_recal_file,tranches_file)
        tool.judge_then_exec(run_sample_id,cmd_recal,snps_recal_file)

        snps_VQSR_file = dir_VQSR + sample_name + ".HC.snps.VQSR.vcf.gz"
        cmd_snps = "{} ApplyVQSR --tmp-dir {} -R {} -variant {} -ts-filter-level 99.5 -tranches-file {} -recal-file {} -mode SNP -O {}"\
        .format(gatk,tmp_dir,ref,vcf_file,tranches_file,snps_recal_file,snps_VQSR_file)
        tool.judge_then_exec(run_sample_id,cmd_snps,snps_VQSR_file)

        ### Indel Recalibration
        indel_recal_file = dir_VQSR + sample_name + ".HC.snps.indels.recal"
        indel_tranches_file = dir_VQSR + sample_name + ".HC.snps.indels.tranches"
        bundle_file1 = bundle_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        indel_recal_cmd = "{} VariantRecalibrator --tmp-dir {} -R {} --variant {} -resource:mills,known=true,training=true,truth=true,prior=12.0 {} -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL --output {} -tranches-file {}"\
        .format(gatk,tmp_dir,ref,snps_VQSR_file,bundle_file1,indel_recal_file,indel_tranches_file)
        tool.judge_then_exec(run_sample_id,indel_recal_cmd,indel_recal_file)

        indels_snps_VQSR_file = dir_VQSR + sample_name + ".HC.snps.indels.VQSR.vcf"
        indels_snps_cmd = "{} ApplyVQSR --tmp-dir {} -R {} -variant {} -ts-filter-level 99.0 -tranches-file {} -recal-file {} -mode INDEL -O {}"\
        .format(gatk,tmp_dir,ref,snps_VQSR_file,indel_tranches_file,indel_recal_file,indels_snps_VQSR_file)
        tool.judge_then_exec(run_sample_id,indels_snps_cmd,indels_snps_VQSR_file)

        indels_snps_VQSR_pass_file = dir_VQSR + sample_name + ".HC.snps.indels.VQSR.PASS.vcf"
        cmd_pass_indel = "awk -F \'\\t\' \'{if($0 ~ /\#/) print; else if($7 == \"PASS\") print}\'" +" {} > {}"\
        .format(indels_snps_VQSR_file,indels_snps_VQSR_pass_file)
        tool.judge_then_exec(run_sample_id,cmd_pass_indel,indels_snps_VQSR_pass_file)

        cmd_gz_1 = "gzip -f {}".format(indels_snps_VQSR_file)
        cmd_gz_2 = "gzip -f {}".format(indels_snps_VQSR_pass_file)
        tool.exec_cmd(cmd_gz_1,run_sample_id)
        tool.exec_cmd(cmd_gz_2,run_sample_id)

    return dir_pre_vcf,dir_VQSR,BQSR_index_rmdup_sorted_bam_filename


def Mutect2(run_sample_id,sample_name,tool,configure,pathes):
    ref = pathes['path']['hg38_ref_file']
    tmp_dir = configure['path']['tmp_dir'] + "/"  # Temporary files directory
    java17 = pathes['path']['java17']
    gatk_jar = pathes['path']['gatk_4.6_jar']
    mem = configure['args']['mem']
    bed_file = configure['others']['bed_file']  # Target regions BED file
    gatk = "{} -Xmx{} -jar {}".format(java17,mem,gatk_jar)

    tumor_with_matched_normal = configure['others']['tumor_with_matched_normal']
    pre_mutation_calling_only = configure['others']['pre_mutation_calling_only']
    if tumor_with_matched_normal:
        tumor_sample = sample_name.split(",")[0]
        normal_sample = sample_name.split(",")[1]
        dir_pre_vcf_tumor,dir_VQSR_tumor,tumor_bam = mutation_calling(run_sample_id,tumor_sample,tool,configure,pathes)
        dir_pre_vcf_normal,dir_VQSR_normal,normal_bam = mutation_calling(run_sample_id,normal_sample,tool,configure,pathes)
        
        if pre_mutation_calling_only:
            return
        
        ### Tumor with matched normal
        mutect_file = dir_pre_vcf_tumor + tumor_sample + '.mutect.vcf.gz'
        cmd_mutect = "{} Mutect2 --tmp-dir {} -R {} -I {} -I {} -normal {} -O {} -L {}"\
        .format(gatk,tmp_dir,ref,tumor_bam,normal_bam,normal_sample,mutect_file,bed_file)
        tool.judge_then_exec(run_sample_id,cmd_mutect,mutect_file)

        filtered_mutect_file = dir_VQSR_tumor + tumor_sample + '.mutect.filtered.vcf.gz'
        cmd_mutect_filter = "{} FilterMutectCalls --tmp-dir {} -R {} -V {} -O {}"\
        .format(gatk,tmp_dir,ref,mutect_file,filtered_mutect_file)
        tool.judge_then_exec(run_sample_id,cmd_mutect_filter,filtered_mutect_file)

    else:
        dir_pre_vcf,dir_VQSR,BQSR_index_rmdup_sorted_bam_filename  = mutation_calling(sample_name,sample_name,tool,configure,pathes)
        
        if pre_mutation_calling_only:
            return
        
        ### Tumor-only mode workflow
        mutect_file = dir_pre_vcf + sample_name + '.mutect.vcf.gz'
        cmd_mutect = "{} Mutect2 --tmp-dir {} -R {} -I {} -O {} -L {}"\
        .format(gatk,tmp_dir,ref,BQSR_index_rmdup_sorted_bam_filename,mutect_file,bed_file)
        tool.judge_then_exec(sample_name,cmd_mutect,mutect_file)

        filtered_mutect_file = dir_VQSR + sample_name + '.mutect.filtered.vcf.gz'
        cmd_mutect_filter = "{} FilterMutectCalls --tmp-dir {} -R {} -V {} -O {}"\
        .format(gatk,tmp_dir,ref,mutect_file,filtered_mutect_file)
        tool.judge_then_exec(sample_name,cmd_mutect_filter,filtered_mutect_file)
    

def wes_human_mutation_calling_start(sample_name,tool,configure,pathes):
    mutation_calling_tool = configure['others']['mutation_calling_tool']
    
    if mutation_calling_tool == 'Mutect2':
        Mutect2(sample_name,sample_name,tool,configure,pathes)
    elif mutation_calling_tool == 'HaplotypeCaller':
        mutation_calling(sample_name,sample_name,tool,configure,pathes)
    else:
        tool.write_log(f"Unsupported mutation calling tool","error")
        return
