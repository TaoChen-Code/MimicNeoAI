# coding=utf-8
import os
from datetime import datetime

from mimicneoai.microbial_pipeline.scripts.microbial_abundance import get_data
from mimicneoai.microbial_pipeline.scripts.get_data_for_binding_pred import get_data_for_binding_pred
from mimicneoai.microbial_pipeline.scripts.hla_binding_pred import pvacbind


def HostSequencesRemoving(sample, configure, paths, tool):
    """Remove host sequences by aligning reads to hg38 and T2T references, then collect unmapped reads.

    Args:
        sample (str): Sample ID.
        configure (dict): Pipeline configuration dictionary.
        paths (dict): Paths dictionary (references, tools).
        tool: Execution tool object providing judge_then_exec/exec_cmd/write_log APIs.
    """
    sample = str(sample)

    # ---- I/O & runtime parameters ----
    input_path = configure['path']['input_dir'] + "/"      # raw/QC input root
    output_path = configure['path']['output_dir'] + "/"    # pipeline output root
    thread = configure['args']['thread']                   # threads per sample
    mem_perthread = configure['args']['mem_perthread']     # samtools sort -m
    pair = configure['others']['pair']                     # paired-end flag
    seq_type = configure['others']['seq_type']             # file name suffix
    QC = configure['others']['QC']                         # whether QC step produced QC fastqs

    # ---- Reference genomes ----
    host_fa_hg38 = paths['database']['microbial']['HOST']['HG38']['FA']
    host_fa_t2t  = paths['database']['microbial']['HOST']['T2T']['FA']

    # ---- Step directories ----
    step_name_qc   = configure['step_name']['QC']
    step_name_hg38 = configure['step_name']['hg38']
    step_name_t2t  = configure['step_name']['t2t']

    output_qc   = output_path + f'{sample}/{step_name_qc}/'
    output_hg38 = output_path + f'{sample}/{step_name_hg38}/'
    output_t2t  = output_path + f'{sample}/{step_name_t2t}/'

    # ---- Intermediate file paths (hg38 stage) ----
    sam_file         = f"{output_hg38}{sample}_{seq_type}_hg38_align.sam"
    bam_file         = f"{output_hg38}{sample}_{seq_type}_hg38_align.bam"
    sorted_bam_file  = f"{output_hg38}{sample}_{seq_type}_hg38_align_sort.bam"
    hg38_unmap_bam   = f"{output_hg38}{sample}_{seq_type}_hg38_unmap.bam"
    hg38_unmap_R1_fq = f"{output_hg38}{sample}_{seq_type}_hg38_unmap.R1.fq"
    hg38_unmap_R2_fq = f"{output_hg38}{sample}_{seq_type}_hg38_unmap.R2.fq"
    hg38_unmap_fq    = f"{output_hg38}{sample}_{seq_type}_hg38_unmap.fq"

    # ---- T2T stage outputs ----
    hg38_unmap_t2t_align_sort_bam = f"{output_t2t}{sample}_{seq_type}_hg38_unmap_t2t_align_sort.bam"
    hg38_unmap_t2t_unmap_bam      = f"{output_t2t}{sample}_{seq_type}_hg38_unmap_t2t_unmap.bam"

    # Ensure hg38 output directory exists
    tool.judge_then_exec(sample, f"mkdir -p {output_hg38}", output_hg38)

    # If final unmapped BAM is missing/empty, run the full removal pipeline
    if not os.path.exists(hg38_unmap_t2t_unmap_bam) or os.path.getsize(hg38_unmap_t2t_unmap_bam) == 0:
        # ---- Align to hg38 ----
        if pair:
            if QC:
                # Use QC fastqs from the QC step directory
                in_root = output_qc
                cmd_align_hg38 = (
                    f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_hg38} "
                    f"{in_root}/{sample}/{sample}.QC.R1.fq.gz {in_root}/{sample}/{sample}.QC.R2.fq.gz > {sam_file}"
                )
            else:
                in_root = input_path
                cmd_align_hg38 = (
                    f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_hg38} "
                    f"{in_root}/{sample}/{sample}.R1.fq.gz {in_root}/{sample}/{sample}.R2.fq.gz > {sam_file}"
                )
        else:
            if QC:
                in_root = output_qc
                cmd_align_hg38 = (
                    f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_hg38} "
                    f"{in_root}/{sample}/{sample}.QC.fq.gz > {sam_file}"
                )
            else:
                in_root = input_path
                cmd_align_hg38 = (
                    f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_hg38} "
                    f"{in_root}/{sample}/{sample}.fq.gz > {sam_file}"
                )

        cmd_sam2bam   = f"samtools view -bS -@ {thread} {sam_file} > {bam_file}"
        cmd_rm_sam    = f"rm {sam_file}"
        cmd_sort_bam  = f"samtools sort -@ {thread} -m {mem_perthread} -o {sorted_bam_file} {bam_file}"
        cmd_rm_bam    = f"rm {bam_file}"
        cmd_flagstat1 = f"samtools flagstat -@ {thread} {sorted_bam_file} > {sorted_bam_file}.flagstat.txt"
        cmd_unmapped  = f"samtools view -b -@ {thread} -f 4 -o {hg38_unmap_bam} {sorted_bam_file}"

        tool.judge_then_exec(sample, cmd_align_hg38, sam_file)
        tool.judge_then_exec(sample, cmd_sam2bam, bam_file)
        if os.path.exists(bam_file) and os.path.getsize(bam_file) > 0:
            tool.exec_cmd(cmd_rm_sam, sample)
        tool.judge_then_exec(sample, cmd_sort_bam, sorted_bam_file)
        if os.path.exists(sorted_bam_file) and os.path.getsize(sorted_bam_file) > 0:
            tool.exec_cmd(cmd_rm_bam, sample)

        tool.judge_then_exec(sample, cmd_flagstat1, f"{sorted_bam_file}.flagstat.txt")
        tool.judge_then_exec(sample, cmd_unmapped, hg38_unmap_bam)

        # Convert unmapped hg38 BAM to FASTQ(s)
        if pair:
            cmd_bam2fq = (
                f"samtools fastq -@ {thread} {hg38_unmap_bam} "
                f"-1 {hg38_unmap_R1_fq} -2 {hg38_unmap_R2_fq} -s /dev/null"
            )
            tool.judge_then_exec(sample, cmd_bam2fq, hg38_unmap_R1_fq)
        else:
            cmd_bam2fq = f"samtools fastq -@ {thread} {hg38_unmap_bam} > {hg38_unmap_fq}"
            tool.judge_then_exec(sample, cmd_bam2fq, hg38_unmap_fq)

        # ---- Align the hg38-unmapped reads to T2T and collect unmapped ----
        tool.judge_then_exec(sample, f"mkdir -p {output_t2t}", output_t2t)

        if pair:
            cmd_t2t_align_sort = (
                f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_t2t} "
                f"{hg38_unmap_R1_fq} {hg38_unmap_R2_fq} | "
                f"samtools view -bS -@ {thread} | "
                f"samtools sort -@ {thread} -m {mem_perthread} -o {hg38_unmap_t2t_align_sort_bam}"
            )
        else:
            cmd_t2t_align_sort = (
                f"bwa mem -q -t {thread} -R '@RG\\tID:{sample}\\tSM:{sample}' {host_fa_t2t} "
                f"{hg38_unmap_fq} | "
                f"samtools view -bS -@ {thread} | "
                f"samtools sort -@ {thread} -m {mem_perthread} -o {hg38_unmap_t2t_align_sort_bam}"
            )
        tool.judge_then_exec(sample, cmd_t2t_align_sort, hg38_unmap_t2t_align_sort_bam)

        cmd_t2t_unmapped = (
            f"samtools view -b -@ {thread} -f 4 -o {hg38_unmap_t2t_unmap_bam} {hg38_unmap_t2t_align_sort_bam}"
        )
        tool.judge_then_exec(sample, cmd_t2t_unmapped, hg38_unmap_t2t_unmap_bam)

    # Always regenerate stats for the final unmapped T2T BAM (useful for downstream depth checks)
    cmd_flagstat2 = (
        f"samtools flagstat -@ {thread} {hg38_unmap_t2t_unmap_bam} > {hg38_unmap_t2t_unmap_bam}.flagstat.txt"
    )
    tool.judge_then_exec(sample, cmd_flagstat2, f"{hg38_unmap_t2t_unmap_bam}.flagstat.txt")


def MicrobialTaxasQuantification(sample, configure, paths, tool):
    """Quantify microbial taxa using GATK PathSeq and prepare inputs for downstream steps.

    Args:
        sample (str): Sample ID.
        configure (dict): Pipeline configuration dictionary.
        paths (dict): Paths dictionary (references, tools).
        tool: Execution tool object.
    """
    sample = str(sample)

    # ---- Configuration ----
    tmp_dir     = configure['path']['tmp_dir'] + "/"
    output_path = configure['path']['output_dir'] + "/"
    thread      = configure['args']['thread']
    mem         = configure['args']['mem']
    pair        = configure['others']['pair']
    seq_type    = configure['others']['seq_type']
    match_length_threshold = configure['others']['match_length_threshold']
    score_threshold        = configure['others']['score_threshold']

    # ---- Reference databases (unified keys) ----
    microbe_dict          = paths['database']['microbial']['MICROBES']['DICT']
    microbe_img           = paths['database']['microbial']['MICROBES']['IMG']
    taxonomy_file         = paths['database']['microbial']['MICROBES']['TAXONOMY_DB']
    microbial_genome_len  = paths['database']['microbial']['MICROBES']['GENOME_LENGTHS']
    catalog_genome_rank   = paths['database']['microbial']['MICROBES']['CATALOG_GENOME']
    tax_id_hierarchy_file = paths['database']['microbial']['MICROBES']['TAXID_HIERARCHY']

    # ---- Tools ----
    gatk_jar = paths['path']['common']['GATK_JAR']

    # ---- Step directories ----
    step_name_pathseq = configure['step_name']['pathseq']
    step_name_nucleic = configure['step_name']['nucleic']
    step_name_hg38    = configure['step_name']['hg38']
    step_name_t2t     = configure['step_name']['t2t']

    output_pathseq = output_path + f'{sample}/{step_name_pathseq}/'
    output_nucleic = output_path + f'{sample}/{step_name_nucleic}/'
    output_hg38    = output_path + f'{sample}/{step_name_hg38}/'
    output_t2t     = output_path + f'{sample}/{step_name_t2t}/'

    # Final BAM of unmapped-to-both (hg38 -> T2T)
    hg38_unmap_t2t_unmap_bam = f"{output_t2t}{sample}_{seq_type}_hg38_unmap_t2t_unmap.bam"

    # Ensure output directories exist
    tool.judge_then_exec(sample, f"mkdir -p {output_pathseq}", output_pathseq)

    # ---- Run PathSeq (Spark local mode) ----
    # Note: host filter images can be injected via paths if desired; here we run minimal config.
    cmd_pathseq = (
        f"java -Xms{mem} -Xmx{mem} -jar {gatk_jar} PathSeqPipelineSpark "
        f"--spark-master local[{thread}] "
        f"--tmp-dir {tmp_dir} "
        f"--input {hg38_unmap_t2t_unmap_bam} "
        f"--min-clipped-read-length 50 "
        f"--microbe-dict {microbe_dict} "
        f"--microbe-bwa-image {microbe_img} "
        f"--taxonomy-file {taxonomy_file} "
        f"--output {output_pathseq}{sample}_{seq_type}_output.pathseq.bam "
        f"--scores-output {output_pathseq}{sample}_{seq_type}_output.pathseq.txt "
        f"--score-metrics {output_pathseq}{sample}_{seq_type}_score.metrics.txt "
        f"--filter-metrics {output_pathseq}{sample}_{seq_type}_filter.metrics.txt "
        f"--divide-by-genome-length true"
    )
    tool.judge_then_exec(sample, cmd_pathseq, f"{output_pathseq}{sample}_{seq_type}_output.pathseq.txt")

    # Prepare nucleic output and convert BAM -> text for downstream parsing
    tool.judge_then_exec(sample, f"mkdir -p {output_nucleic}", output_nucleic)

    cmd_bam2txt = (
        f"samtools view -@ {thread} {output_pathseq}{sample}_{seq_type}_output.pathseq.bam "
        f"> {output_nucleic}{sample}.txt"
    )
    tool.judge_then_exec(sample, cmd_bam2txt, f"{output_nucleic}{sample}.txt")

    # ---- Extract microbial sequences for BLASTX & abundance ----
    # Skip heavy recomputation if final genus CSV already exists
    if not os.path.exists(f"{output_nucleic}/{sample}.microbe_abundance_genus.csv"):
        tool.write_log(f"[{sample}] Extracting microbial sequences for BLASTX", "info")
        start = datetime.now()
        get_data(
            sample,
            catalog_genome_rank_file=catalog_genome_rank,
            tax_id_hierarchy_file=tax_id_hierarchy_file,
            output_hg38=output_hg38,
            output_nucleic=output_nucleic,
            microbial_genome_length=microbial_genome_len,
            match_length_threshold=match_length_threshold,
            score_threshold=score_threshold,
            pair=pair,
            seq_type=seq_type,
        )
        end = datetime.now()
        using_time = tool.print_time(end - start)
        tool.write_log(f"[{sample}] Sequence extraction completed in {using_time}", "info")


def MicrobialPeptidesIdentification(sample, configure, paths, tool):
    """Identify microbial peptides via BLASTX against a protein database.

    Args:
        sample (str): Sample ID.
        configure (dict): Pipeline configuration dictionary.
        paths (dict): Paths dictionary (references, tools).
        tool: Execution tool object.
    """
    output_path = configure['path']['output_dir'] + "/"
    step_name_nucleic = configure['step_name']['nucleic']
    step_name_blastx  = configure['step_name']['blastx']

    output_nucleic = output_path + f'{sample}/{step_name_nucleic}/'
    output_blastx  = output_path + f'{sample}/{step_name_blastx}/'

    thread = configure['args']['thread']
    min_pident_length = configure['others']['min_pident_length']

    # BLAST database config
    db_dir = paths['database']['microbial']['BLAST']['DB_DIR']
    outfmt = paths['database']['microbial']['BLAST']['OUTFMT']

    # Ensure output directory exists
    tool.judge_then_exec(sample, f"mkdir -p {output_blastx}", output_blastx)

    # Run BLASTX on nucleic sequences produced earlier
    blastx_out = f"{output_blastx}{sample}.taxidlist.blast"
    cmd_blastx = (
        f"blastx -query {output_nucleic}{sample}.fasta "
        f"-out {blastx_out} "
        f"-db {db_dir} -outfmt {outfmt} -evalue 5e-2 "
        f"-num_threads {thread} -max_target_seqs 5"
    )
    tool.judge_then_exec(sample, cmd_blastx, blastx_out)

    # Convert BLASTX hits to pVACbind-ready peptide FASTA if not already created
    peptide_fa = f"{output_blastx}/{sample}.peptide.fasta"
    if not os.path.exists(peptide_fa):
        tool.write_log("[INFO] Processing BLASTX results for pVACbind input", "info")
        start = datetime.now()
        get_data_for_binding_pred(
            f"{sample}.taxidlist.blast",
            f"{sample}.peptide.fasta",
            output_blastx,
            min_pident_length
        )
        end = datetime.now()
        tool.print_time(end - start)


def MicrobialPeptidesBindingPrediction(sample, configure, paths, tool):
    """Run pVACbind for peptide–MHC binding prediction on BLASTX-derived peptides.

    Args:
        sample (str): Sample ID.
        configure (dict): Pipeline configuration dictionary.
        paths (dict): Paths dictionary (references, tools).
        tool: Execution tool object.
    """
    output_path = configure['path']['output_dir'] + "/"
    step_name_blastx   = configure['step_name']['blastx']
    step_name_hla      = configure['step_name']['hla']
    step_name_pvacbind = configure['step_name']['pvacbind']

    output_blastx   = output_path + f'{sample}/{step_name_blastx}/'
    output_pvacbind = output_path + f'{sample}/{step_name_pvacbind}/'

    # Only run if peptide FASTA exists
    peptide_fa = f"{output_blastx}/{sample}.peptide.fasta"
    if os.path.exists(peptide_fa):
        tool.judge_then_exec(sample, f"mkdir -p {output_pvacbind}", output_pvacbind)
        pvacbind(sample, configure, paths, tool)
