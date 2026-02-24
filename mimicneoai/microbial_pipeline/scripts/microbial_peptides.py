# coding=utf-8
import os
from datetime import datetime
from mimicneoai.functions.utils import format_java_heap
import pandas as pd
from mimicneoai.microbial_pipeline.scripts.get_data_for_blastx import get_data
from mimicneoai.microbial_pipeline.scripts.get_data_for_binding_pred import get_data_for_binding_pred
from mimicneoai.microbial_pipeline.scripts.hla_binding_pred import pvacbind

def HostSequencesRemoving(sample, configure, paths, tool):
    """
    Remove host sequences by aligning reads to hg38 then T2T, and collecting unmapped reads.

    Pipeline (paired-end, strict):
      FASTQ -> bwa mem hg38 -> name-sorted BAM -> extract (-f 12) -> FASTQ (R1/R2 + single)
            -> bwa mem T2T  -> name-sorted BAM -> extract (-f 12) -> final unmapped BAM
      Always: samtools flagstat on final unmapped BAM.
    """
    sample = str(sample)

    # ------------------------------------------------------------------
    # 0) Read configs / references
    # ------------------------------------------------------------------
    input_root  = configure["path"]["input_dir"].rstrip("/")   # raw/QC input root
    output_root = configure["path"]["output_dir"].rstrip("/")  # pipeline output root

    thread       = int(configure["args"]["thread"])
    mem_perthread = configure["args"]["mem_perthread"]

    pair     = bool(configure["others"]["pair"])
    seq_type = str(configure["others"]["seq_type"])
    QC       = bool(configure["others"]["QC"])

    host_fa_hg38 = paths["database"]["microbial"]["HOST"]["HG38"]["FA"]
    host_fa_t2t  = paths["database"]["microbial"]["HOST"]["T2T"]["FA"]

    step_qc   = configure["step_name"]["QC"]
    step_hg38 = configure["step_name"]["hg38"]
    step_t2t  = configure["step_name"]["t2t"]

    # ------------------------------------------------------------------
    # 1) Utilities (centralize command construction)
    # ------------------------------------------------------------------
    def _mkdir(p: str):
        tool.judge_then_exec(sample, f"mkdir -p {p}", p)

    def _rg() -> str:
        # Note: keep exact escaping required by bwa
        return f"@RG\\tID:{sample}\\tSM:{sample}"

    # def _unmapped_flag(is_paired: bool) -> int:
    #     # paired: 0x4 (read unmapped) + 0x8 (mate unmapped) = 12, strict both ends unmapped
    #     return 12 if is_paired else 4

    def _unmapped_flag() -> int:
        # Always keep reads where the read itself is unmapped (0x4)
        return 4

    def _flagstat_cmd(bam: str) -> str:
        return f"samtools flagstat -@ {thread} {bam} > {bam}.flagstat.txt"

    def _fastq_inputs():
        """
        Return FASTQ input paths depending on QC/pair setting.
        Keeps your existing directory convention (note the extra {sample}/).
        """
        if QC:
            qc_dir = f"{output_root}/{sample}/{step_qc}".rstrip("/")
            # Your original code uses: {output_qc}/{sample}/{sample}.QC.R1.fq.gz
            base = f"{qc_dir}/{sample}/{sample}.QC"
        else:
            base = f"{input_root}/{sample}/{sample}"

        if pair:
            return (f"{base}.R1.fq.gz", f"{base}.R2.fq.gz")
        else:
            return (f"{base}.fq.gz",)

    # ------------------------------------------------------------------
    # 2) Paths (grouped, predictable)
    # ------------------------------------------------------------------
    out_hg38_dir = f"{output_root}/{sample}/{step_hg38}/"
    out_t2t_dir  = f"{output_root}/{sample}/{step_t2t}/"

    # hg38 intermediate
    hg38_name_bam   = f"{out_hg38_dir}{sample}_{seq_type}.hg38.name_sorted.bam"
    hg38_unmap_bam  = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.bam"

    # hg38 unmapped fastqs
    hg38_unmap_r1   = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.R1.fq"
    hg38_unmap_r2   = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.R2.fq"
    hg38_unmap_anom    = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.anomalous.fq"   # -0
    hg38_unmap_single  = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.singleton.fq"   # -s (NEW)
    hg38_unmap_se   = f"{out_hg38_dir}{sample}_{seq_type}.hg38.unmapped.fq"

    # T2T outputs
    t2t_name_bam    = f"{out_t2t_dir}{sample}_{seq_type}.hg38unmap.t2t.name_sorted.bam"
    final_unmap_bam = f"{out_t2t_dir}{sample}_{seq_type}.hg38unmap.t2t.unmapped.bam"

    _mkdir(out_hg38_dir)
    _mkdir(out_t2t_dir)

    # ------------------------------------------------------------------
    # 3) Decide whether to run full pipeline
    # ------------------------------------------------------------------
    need_run = (not os.path.exists(final_unmap_bam)) or (os.path.getsize(final_unmap_bam) == 0)

    # ------------------------------------------------------------------
    # 4) hg38 stage: FASTQ -> hg38 align -> name-sort -> extract unmapped -> FASTQ
    # ------------------------------------------------------------------
    if need_run:
        fastqs_in = _fastq_inputs()

        # 4.1–4.3: name-sorted BAM
        fq_part = " ".join(fastqs_in)
        cmd_align_sort_hg38 = (
            f"bwa mem -q -t {thread} "
            f"-R '{_rg()}' "
            f"{host_fa_hg38} {fq_part} | "
            f"samtools view -b -@ {thread} - | "
            f"samtools sort -n -@ {thread} -m {mem_perthread} "
            f"-o {hg38_name_bam} -"
        )
        tool.judge_then_exec(sample, cmd_align_sort_hg38, hg38_name_bam)

        # 4.4 stats on name-sorted BAM
        tool.judge_then_exec(sample, _flagstat_cmd(hg38_name_bam), f"{hg38_name_bam}.flagstat.txt")


        # 4.5 extract unmapped
        flag = _unmapped_flag()
        cmd_extract_hg38_unmap = (
            f"samtools view -b -@ {thread} -f {flag} "
            f"-o {hg38_unmap_bam} {hg38_name_bam}"
        )
        tool.judge_then_exec(sample, cmd_extract_hg38_unmap, hg38_unmap_bam)

        # 4.6 BAM -> FASTQ(s) for T2T
        if pair:
            # -1/-2: proper pairs
            # -0: anomalous/other reads (retained)
            # -s: singleton reads (NEW: retained, previously /dev/null)
            cmd_bam2fq = (
                f"samtools fastq -@ {thread} {hg38_unmap_bam} "
                f"-1 {hg38_unmap_r1} "
                f"-2 {hg38_unmap_r2} "
                f"-0 {hg38_unmap_anom} "
                f"-s {hg38_unmap_single}"
            )
            tool.judge_then_exec(sample, cmd_bam2fq, hg38_unmap_r1)
        else:
            cmd_bam2fq = f"samtools fastq -@ {thread} {hg38_unmap_bam} > {hg38_unmap_se}"
            tool.judge_then_exec(sample, cmd_bam2fq, hg38_unmap_se)

        # ------------------------------------------------------------------
        # 5) T2T stage: hg38-unmapped FASTQ -> T2T align+name-sort -> extract unmapped
        # ------------------------------------------------------------------
        if pair:
            t2t_fastqs = (hg38_unmap_r1, hg38_unmap_r2)
        else:
            t2t_fastqs = (hg38_unmap_se,)

        fq_t2t_part = " ".join(t2t_fastqs)
        cmd_align_sort_t2t = (
            f"bwa mem -q -t {thread} "
            f"-R '{_rg()}' "
            f"{host_fa_t2t} {fq_t2t_part} | "
            f"samtools view -b -@ {thread} - | "
            f"samtools sort -n -@ {thread} -m {mem_perthread} "
            f"-o {t2t_name_bam} -"
        )
        tool.judge_then_exec(sample, cmd_align_sort_t2t, t2t_name_bam)


        t2t_flag = _unmapped_flag()
        cmd_extract_t2t_unmap = (
            f"samtools view -b -@ {thread} -f {t2t_flag} "
            f"-o {final_unmap_bam} {t2t_name_bam}"
        )
        tool.judge_then_exec(sample, cmd_extract_t2t_unmap, final_unmap_bam)

    # ------------------------------------------------------------------
    # 6) Always: stats for final output
    # ------------------------------------------------------------------
    tool.judge_then_exec(sample, _flagstat_cmd(final_unmap_bam), f"{final_unmap_bam}.flagstat.txt")

def VectorContaminationRemoving(
    sample: str,
    configure: dict,
    paths: dict,
    tool,
):
    """
    Remove vector contamination (UniVec) from host-removed BAM.

    Input BAM (fixed):
      <output_dir>/<sample>/<step_t2t>/<sample>_<seq_type>.hg38unmap.t2t.unmapped.bam

    Output:
      <output_dir>/<sample>/<step_vector>/
        ├── 01.fastq/
        │     ├── *.R1.fq.gz
        │     ├── *.R2.fq.gz
        │     ├── *.singleton.fq.gz
        │     └── *.cat0.fq.gz
        └── 02.bam/
              ├── *.paired.vector_unmapped.bam
              ├── *.single.vector_unmapped.bam
              └── *.vector_unmapped.merged.bam   (FINAL)
    """
    sample = str(sample)

    # ------------------------------------------------------------------
    # 0) Read configs
    # ------------------------------------------------------------------
    output_root = configure["path"]["output_dir"].rstrip("/")
    thread = int(configure["args"]["thread"])
    mem_perthread = configure["args"].get("mem_perthread", "1G")

    seq_type = str(configure["others"]["seq_type"])

    step_t2t    = configure["step_name"]["t2t"]
    step_vector = configure["step_name"]["vector"]

    # ------------------------------------------------------------------
    # 1) Input BAM (fixed name from HostSequencesRemoving)
    # ------------------------------------------------------------------
    in_bam = (
        f"{output_root}/{sample}/{step_t2t}/"
        f"{sample}_{seq_type}.hg38unmap.t2t.unmapped.bam"
    )

    # ------------------------------------------------------------------
    # 2) UniVec reference
    # ------------------------------------------------------------------
    univec_ref = paths["database"]["microbial"]["MICROBES"]["UniVec"]
    if os.path.isdir(univec_ref):
        vector_fa = os.path.join(univec_ref, "UniVec.fa")
    else:
        vector_fa = univec_ref

    # ------------------------------------------------------------------
    # 3) Output directories
    # ------------------------------------------------------------------
    out_dir = f"{output_root}/{sample}/{step_vector}/"
    fq_dir  = f"{out_dir}01.fastq/"
    bam_dir = f"{out_dir}02.bam/"

    tool.judge_then_exec(sample, f"mkdir -p {fq_dir}", fq_dir)
    tool.judge_then_exec(sample, f"mkdir -p {bam_dir}", bam_dir)

    # ------------------------------------------------------------------
    # 4) Filenames
    # ------------------------------------------------------------------
    r1_fq = f"{fq_dir}{sample}_{seq_type}.host_removed.R1.fq.gz"
    r2_fq = f"{fq_dir}{sample}_{seq_type}.host_removed.R2.fq.gz"
    se_fq = f"{fq_dir}{sample}_{seq_type}.host_removed.singleton.fq.gz"
    cat0_fq = f"{fq_dir}{sample}_{seq_type}.host_removed.cat0.fq.gz"

    pe_unmap_bam = f"{bam_dir}{sample}_{seq_type}.paired.vector_unmapped.bam"
    se_unmap_bam = f"{bam_dir}{sample}_{seq_type}.single.vector_unmapped.bam"
    merged_bam   = f"{bam_dir}{sample}_{seq_type}.vector_unmapped.merged.bam"

    # ------------------------------------------------------------------
    # 5) BAM → FASTQ (PE / singleton / cat0)
    # ------------------------------------------------------------------
    cmd_bam2fq = (
        f"samtools sort -n -@ {thread} -m {mem_perthread} {in_bam} | "
        f"samtools fastq -@ {thread} - "
        f"-1 {r1_fq} -2 {r2_fq} "
        f"-s {se_fq} -0 {cat0_fq}"
    )
    tool.judge_then_exec(sample, cmd_bam2fq, r1_fq)

    # ------------------------------------------------------------------
    # 6) BWA → UniVec (PE)
    # ------------------------------------------------------------------
    pe_vec_bam = f"{bam_dir}{sample}_{seq_type}.paired.vector.bam"
    cmd_bwa_pe = (
        f"bwa mem -q -t {thread} "
        f"-R '@RG\\tID:{sample}\\tSM:{sample}' "
        f"{vector_fa} {r1_fq} {r2_fq} | "
        f"samtools view -b -@ {thread} -o {pe_vec_bam} -"
    )
    tool.judge_then_exec(sample, cmd_bwa_pe, pe_vec_bam)

    # ------------------------------------------------------------------
    # 7) BWA → UniVec (SE)
    # ------------------------------------------------------------------
    se_vec_bam = f"{bam_dir}{sample}_{seq_type}.single.vector.bam"
    cmd_bwa_se = (
        f"bwa mem -q -t {thread} "
        f"-R '@RG\\tID:{sample}_SE\\tSM:{sample}' "
        f"{vector_fa} {se_fq} | "
        f"samtools view -b -@ {thread} -o {se_vec_bam} -"
    )
    tool.judge_then_exec(sample, cmd_bwa_se, se_vec_bam)

    # ------------------------------------------------------------------
    # 8) Extract vector-unmapped reads (统一 -f 4)
    # ------------------------------------------------------------------
    unmapped_flag = 4

    cmd_unmap_pe = (
        f"samtools view -b -@ {thread} -f {unmapped_flag} "
        f"-o {pe_unmap_bam} {pe_vec_bam}"
    )
    cmd_unmap_se = (
        f"samtools view -b -@ {thread} -f {unmapped_flag} "
        f"-o {se_unmap_bam} {se_vec_bam}"
    )

    tool.judge_then_exec(sample, cmd_unmap_pe, pe_unmap_bam)
    tool.judge_then_exec(sample, cmd_unmap_se, se_unmap_bam)

    # ------------------------------------------------------------------
    # 9) Merge PE / SE unmapped BAMs
    # ------------------------------------------------------------------
    have_pe = os.path.exists(pe_unmap_bam) and os.path.getsize(pe_unmap_bam) > 0
    have_se = os.path.exists(se_unmap_bam) and os.path.getsize(se_unmap_bam) > 0

    if have_pe and have_se:
        cmd_merge = f"samtools merge -@ {thread} {merged_bam} {pe_unmap_bam} {se_unmap_bam}"
        tool.judge_then_exec(sample, cmd_merge, merged_bam)
    elif have_pe:
        cmd_copy = f"samtools view -b {pe_unmap_bam} -o {merged_bam}"
        tool.judge_then_exec(sample, cmd_copy, merged_bam)
    elif have_se:
        cmd_copy = f"samtools view -b {se_unmap_bam} -o {merged_bam}"
        tool.judge_then_exec(sample, cmd_copy, merged_bam)
    else:
        tool.write_log(sample, f"[WARN] {sample}: no vector-unmapped reads.\n")

    # ------------------------------------------------------------------
    # 10) QC stats (optional but recommended)
    # ------------------------------------------------------------------
    if os.path.exists(merged_bam) and os.path.getsize(merged_bam) > 0:
        cmd_flagstat = (
            f"samtools flagstat -@ {thread} {merged_bam} "
            f"> {merged_bam}.flagstat.txt"
        )
        tool.judge_then_exec(sample, cmd_flagstat, f"{merged_bam}.flagstat.txt")



def MicrobialTaxasQuantification(sample, configure, paths, tool):
    """Quantify microbial taxa using GATK PathSeq and prepare inputs for downstream steps.

    Notes (updated):
      - Host removal (hg38 + T2T) has already been done upstream.
      - Vector contamination removal (UniVec) has already been done upstream.
      - Therefore PathSeq is run WITHOUT host filtering (--filter-bwa-image / --kmer-file removed),
        and takes the vector-clean merged BAM as input.
    """
    sample = str(sample)

    # ---- Configuration ----
    tmp_dir     = configure['path']['tmp_dir'].rstrip("/") + "/"
    output_path = configure['path']['output_dir'].rstrip("/") + "/"
    thread      = int(configure['args']['thread'])
    # PathSeqPipelineSpark 建议限制本地线程数，避免资源过度竞争
    pathseq_threads = min(thread, 5)

    mem = configure['args']['mem']
    # ---- Example usage ----
    _, xmx, xms = format_java_heap(configure["args"]["mem"], xms_ratio=0.5)

    pair     = configure['others']['pair']
    seq_type = configure['others']['seq_type']

    # ---- Reference databases ----
    microbe_dict  = paths['database']['microbial']['MICROBES']['DICT']
    microbe_img   = paths['database']['microbial']['MICROBES']['IMG']
    taxonomy_file = paths['database']['microbial']['MICROBES']['TAXONOMY_DB']

    # ---- Tools ----
    gatk_jar = paths['path']['common']['GATK_JAR']

    # ---- Step directories ----
    step_name_pathseq = configure['step_name']['pathseq']
    step_name_nucleic = configure['step_name']['nucleic']
    step_name_vector  = configure['step_name']['vector']  # NEW: vector step

    output_pathseq = output_path + f'{sample}/{step_name_pathseq}/'
    output_nucleic = output_path + f'{sample}/{step_name_nucleic}/'

    # Vector step BAM dir (matches your VectorContaminationRemoving layout)
    output_vector = output_path + f'{sample}/{step_name_vector}/'
    bam_dir = output_vector + "02.bam/"

    # ---- Input BAM: vector-unmapped merged BAM ----
    merged_bam = f"{bam_dir}{sample}_{seq_type}.vector_unmapped.merged.bam"

    # Ensure output directories exist
    tool.judge_then_exec(sample, f"mkdir -p {output_pathseq}", output_pathseq)
    tool.judge_then_exec(sample, f"mkdir -p {output_nucleic}", output_nucleic)

    # Validate input BAM exists
    if not os.path.exists(merged_bam) or os.path.getsize(merged_bam) == 0:
        raise FileNotFoundError(
            f"[{sample}] Vector-clean merged BAM is missing/empty: {merged_bam}"
        )

    pathseq_bam     = f"{output_pathseq}{sample}_{seq_type}_output.pathseq.bam"
    scores_txt      = f"{output_pathseq}{sample}_{seq_type}_output.pathseq.txt"
    score_metrics   = f"{output_pathseq}{sample}_{seq_type}_score.metrics.txt"
    filter_metrics  = f"{output_pathseq}{sample}_{seq_type}_filter.metrics.txt"

    # ---- Run PathSeq (Spark local mode) ----
    # If PathSeq scores already exist, skip rerun to save time.
    if not os.path.exists(scores_txt) or os.path.getsize(scores_txt) == 0:
        cmd_pathseq = (
            f"java -Xms{xms} -Xmx{xmx} -jar {gatk_jar} PathSeqPipelineSpark "
            f"--spark-master local[{pathseq_threads}] "
            f"--tmp-dir {tmp_dir} "
            f"--input {merged_bam} "
            f"--min-clipped-read-length 50 "
            # REMOVED (host filtering already done upstream):
            # f"--filter-bwa-image {host_img} "
            # f"--kmer-file {host_kmer} "
            f"--microbe-dict {microbe_dict} "
            f"--microbe-bwa-image {microbe_img} "
            f"--taxonomy-file {taxonomy_file} "
            f"--output {pathseq_bam} "
            f"--scores-output {scores_txt} "
            f"--score-metrics {score_metrics} "
            f"--filter-metrics {filter_metrics} "
            f"--divide-by-genome-length true"
        )
        tool.judge_then_exec(sample, cmd_pathseq, scores_txt)

    # === New behavior (unchanged) ===
    # Directly read PathSeq BAM and extract YP-tagged reads into a fasta.gz file,
    # plus a small summary TSV describing read counts.
    summary_tsv = os.path.join(output_nucleic, f"{sample}.pathseq_selected.summary.tsv")

    # Skip re-extraction if summary already exists (idempotent behavior).
    if not os.path.exists(summary_tsv) or os.path.getsize(summary_tsv) == 0:
        tool.write_log(f"[{sample}] Extracting YP-tagged reads from PathSeq BAM", "info")
        start = datetime.now()
        get_data(
            sample,
            output_pathseq=output_pathseq,
            output_nucleic=output_nucleic,
            pair=pair,
            seq_type=seq_type,
        )
        end = datetime.now()
        using_time = tool.print_time(end - start)
        tool.write_log(f"[{sample}] YP-tagged read extraction completed in {using_time}", "info")



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


    # BLAST database config
    db_dir = paths['database']['microbial']['BLAST']['DB_DIR']
    outfmt = paths['database']['microbial']['BLAST']['OUTFMT']

    # Ensure output directory exists
    tool.judge_then_exec(sample, f"mkdir -p {output_blastx}", output_blastx)

    # Run BLASTX on nucleic sequences produced earlier
    fa_gz = f"{output_nucleic}{sample}.pathseq_selected.fa.gz"
    fa_fa = f"{output_nucleic}{sample}.pathseq_selected.fa"
    blastx_out = f"{output_blastx}{sample}.blastx"

    # 1) Decompress the gzipped FASTA.
    #    This step uses exec_cmd instead of judge_then_exec because it is considered
    #    an internal preprocessing step and should not appear as the main “Running:” command.
    cmd_unzip = f"zcat {fa_gz} > {fa_fa}"
    tool.exec_cmd(cmd_unzip, sample)

    # 2) Run BLASTX using the uncompressed FASTA.
    #    judge_then_exec is used here so that “blastx” appears in the main log entry.
    cmd_blastx = (
        f"blastx "
        f"-query {fa_fa} "
        f"-out {blastx_out} "
        f"-db {db_dir} "
        f"-outfmt {outfmt} "
        f"-num_threads {thread} "
    )
    tool.judge_then_exec(sample, cmd_blastx, blastx_out)

    # 3) Optionally remove the temporary uncompressed FASTA to save disk space.
    #    The removal is performed safely via Python rather than `rm -f`,
    #    with additional checks to avoid accidental deletion of unexpected files.
    if os.path.exists(blastx_out) and os.path.getsize(blastx_out) > 0:
        if os.path.isfile(fa_fa):
            # Sanity check to ensure this is the expected temporary FASTA file.
            if fa_fa.endswith(".fa"):
                try:
                    os.remove(fa_fa)
                    tool.write_log(f"[{sample}] Removed temporary FASTA: {fa_fa}", "info")
                except OSError as e:
                    tool.write_log(
                        f"[{sample}] Failed to remove temporary FASTA {fa_fa}: {e}",
                        "warning",
                    )
            else:
                tool.write_log(
                    f"[{sample}] Skipped deletion of {fa_fa}: unexpected extension (not .fa)",
                    "warning",
                )

    # Load catalog table for validating BLASTX hits
    catalog_path = paths['database']['microbial']['BLAST']['CATALOG_PROT']
    catalog_df = pd.read_csv(catalog_path, sep="\t", header=None)
    catalog_df.columns = ['prot_id', 'tax_id']
    catalog_df["tax_id"] = catalog_df["tax_id"].astype(str)

    # Parse BLASTX column names from the configured outfmt string
    # The outfmt string has the form: '6 qseqid qlen sseqid ... qcovhsp qcovs'
    outfmt_cfg = paths['database']['microbial']['BLAST']['OUTFMT']
    outfmt_clean = outfmt_cfg.strip().strip("'\"")  # remove outer quotes
    tokens = outfmt_clean.split()
    if tokens[0] != "6":
        raise ValueError(f"BLAST OUTFMT must start with 6, got: {outfmt_clean}")
    blast_colnames = tokens[1:]  # all field names after the leading "6"

    # Load BLASTX filtering thresholds from the configuration file
    blastx_min_pident = float(configure['others']['blastx_min_percent_identity'])
    blastx_max_evalue = float(configure['others']['blastx_max_evalue'])
    # use BLASTX-reported qcovs as the coverage filter
    blastx_min_qcovs = float(configure['others']['blastx_min_query_coverage'])

    # Convert BLASTX results into a pVACbind-ready FASTA file
    tool.write_log("[INFO] Processing BLASTX results for pVACbind input", "info")
    start = datetime.now()

    get_data_for_binding_pred(
        blast_file=f"{sample}.blastx",
        colnames=blast_colnames,
        pvacbind_file=f"{sample}.peptide.fasta",
        output_blastx=output_blastx,
        min_pident=blastx_min_pident,
        catalog_df=catalog_df,
        max_evalue=blastx_max_evalue,
        min_qcovs=blastx_min_qcovs,   # ← 新参数名
        sample=sample,
        tool=tool,
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
    step_name_pvacbind = configure['step_name']['pvacbind']

    output_blastx   = output_path + f'{sample}/{step_name_blastx}/'
    output_pvacbind = output_path + f'{sample}/{step_name_pvacbind}/'

    # Only run if peptide FASTA exists
    peptide_fa = f"{output_blastx}/{sample}.peptide.fasta"
    if os.path.exists(peptide_fa):
        tool.judge_then_exec(sample, f"mkdir -p {output_pvacbind}", output_pvacbind)
        pvacbind(sample, configure, paths, tool)
