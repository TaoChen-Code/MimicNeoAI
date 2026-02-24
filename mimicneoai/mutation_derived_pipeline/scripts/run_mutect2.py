import os
from typing import Dict, Any
from mimicneoai.functions.utils import format_java_heap


def Mutect2(
    run_sample_id: str,
    sample_name: str,
    tool,
    configure: dict,
    paths: dict,
    *,
    tumor_sample: str,
    normal_sample: str,
    dir_varcall_root_tumor: str,
    tumor_bam: str,
    normal_bam: str,
):
    """
    Matched-normal Mutect2 pipeline:
      Mutect2 -> GetPileupSummaries -> CalculateContamination ->
      LearnReadOrientationModel -> FilterMutectCalls

    Outputs under:
      <output>/<tumor>/04.variants_calling/Mutect2/
    """
    hg38 = paths["database"]["neoantigen"]["HG38"]
    ref = hg38["REF_FASTA"]
    tmp_dir = configure["path"]["tmp_dir"].rstrip("/") + "/"
    gatk_jar = paths["path"]["common"]["GATK_JAR"]
    bed_file = configure["others"]["bed_file"]
    _, xmx, xms = format_java_heap(configure["args"]["mem"], xms_ratio=0.5)
    gatk = f"java -Xms{xms} -Xmx{xmx} -jar {gatk_jar}"

    # Mutect2 resources (uppercase keys in paths.yaml)
    germline_resource = hg38["GERMLINE_RESOURCE_VCF"]
    panel_of_normals = hg38["PANEL_OF_NORMALS_VCF"]
    common_sites_vcf = hg38["COMMON_SITES_VCF"]

    # knobs
    pad_bp = int(configure["others"].get("pad_bp", 100))
    min_bq = int(configure["others"].get("min_base_quality", 20))
    min_af = float(configure["others"].get("min_allele_fraction", 0.02))

    # caller subdir
    dir_mutect2 = os.path.join(dir_varcall_root_tumor, "Mutect2") + "/"
    tool.judge_then_exec(run_sample_id, f"mkdir -p {dir_mutect2}", dir_mutect2)

    mutect_vcf = os.path.join(dir_mutect2, f"{tumor_sample}.mutect.vcf.gz")
    f1r2_tar = mutect_vcf + ".f1r2.tar.gz"
    tumor_pu_table = os.path.join(dir_mutect2, f"{tumor_sample}.tumor.pileups.table")
    normal_pu_table = os.path.join(dir_mutect2, f"{tumor_sample}.normal.pileups.table")
    contam_table = os.path.join(dir_mutect2, f"{tumor_sample}.contamination.table")
    contam_segments = os.path.join(dir_mutect2, f"{tumor_sample}.segments.table")
    ob_model = os.path.join(dir_mutect2, f"{tumor_sample}.read-orientation-model.tar.gz")
    filtered_vcf = os.path.join(dir_mutect2, f"{tumor_sample}.mutect.filtered.vcf.gz")

    # Mutect2
    cmd_mutect = (
        f"{gatk} Mutect2 --tmp-dir {tmp_dir} "
        f"-R {ref} -I {tumor_bam} -I {normal_bam} -normal {normal_sample} "
        f"-O {mutect_vcf} -L {bed_file} --interval-padding {pad_bp} "
        f"--germline-resource {germline_resource} "
        f"--panel-of-normals {panel_of_normals} "
        f"--min-base-quality-score {min_bq} "
        f"--f1r2-tar-gz {f1r2_tar}"
    )
    tool.judge_then_exec(run_sample_id, cmd_mutect, mutect_vcf)

    # contamination
    cmd_pu_tumor = (
        f"{gatk} GetPileupSummaries "
        f"-I {tumor_bam} -V {common_sites_vcf} "
        f"-L {bed_file} --interval-padding {pad_bp} "
        f"-O {tumor_pu_table}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pu_tumor, tumor_pu_table)

    cmd_pu_normal = (
        f"{gatk} GetPileupSummaries "
        f"-I {normal_bam} -V {common_sites_vcf} "
        f"-L {bed_file} --interval-padding {pad_bp} "
        f"-O {normal_pu_table}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pu_normal, normal_pu_table)

    cmd_contam = (
        f"{gatk} CalculateContamination "
        f"-I {tumor_pu_table} -matched {normal_pu_table} "
        f"-O {contam_table} --tumor-segmentation {contam_segments}"
    )
    tool.judge_then_exec(run_sample_id, cmd_contam, contam_table)

    # orientation bias priors
    cmd_ob = f"{gatk} LearnReadOrientationModel -I {f1r2_tar} -O {ob_model}"
    tool.judge_then_exec(run_sample_id, cmd_ob, ob_model)

    # FilterMutectCalls
    cmd_filter = (
        f"{gatk} FilterMutectCalls --tmp-dir {tmp_dir} "
        f"-R {ref} -V {mutect_vcf} -O {filtered_vcf} "
        f"--contamination-table {contam_table} "
        f"--ob-priors {ob_model} "
        f"--min-allele-fraction {min_af}"
    )
    tool.judge_then_exec(run_sample_id, cmd_filter, filtered_vcf)

    # ---- PASS only + index (Mutect2) ----
    # NOTE: Ensure bcftools can be overridden in configure, fallback to PATH
    bcftools = str(configure["others"].get("bcftools", "bcftools")).strip()

    filtered_pass_vcf = os.path.join(dir_mutect2, f"{tumor_sample}.mutect.filtered.PASS.vcf.gz")

    tool.write_log("Mutect2: bcftools PASS + index", "info")
    cmd_pass = (
        f"{bcftools} view -f PASS -Oz -o {filtered_pass_vcf} {filtered_vcf} "
        f"&& {bcftools} index -t {filtered_pass_vcf}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pass, filtered_pass_vcf)

    tool.write_log(
        f"Mutect2 done: {filtered_vcf}; PASS: {filtered_pass_vcf}",
        "info",
    )