import argparse
import os
import pathlib
from typing import Optional


def strelka2_main(argv=None, tool=None, run_sample_id: Optional[str] = None):
    """
    Manta (somatic) -> Strelka2 -> merge -> PASS  (NO VEP)
    All commands executed via tool.judge_then_exec (tool handles logs).
    """
    if tool is None:
        raise ValueError("strelka2_main requires tool")

    ap = argparse.ArgumentParser(
        description="Manta (somatic) -> Strelka2 -> merge -> PASS  (NO VEP)"
    )
    ap.add_argument("--sample-t", required=True, help="Tumor sample name")
    ap.add_argument("--sample-n", required=True, help="Normal sample name")
    ap.add_argument("--ref", required=True, help="Reference FASTA (with .fai)")
    ap.add_argument("--call-regions", required=True, help="Calling BED (.bed or .bed.gz)")
    ap.add_argument("--tumor-bam", required=True, help="Tumor BAM (BQSR)")
    ap.add_argument("--normal-bam", required=True, help="Normal BAM (BQSR)")
    ap.add_argument("--threads", type=int, default=20, help="Threads (default: 20)")
    ap.add_argument("--work-root", required=True, help="Root working dir for outputs")

    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--config-manta", default="configManta.py")
    ap.add_argument("--strelka-config", default="configureStrelkaSomaticWorkflow.py")

    args = ap.parse_args(argv)

    if not run_sample_id:
        run_sample_id = f"{args.sample_t},{args.sample_n}"

    sample_t = args.sample_t
    work_root = pathlib.Path(args.work_root)

    manta_dir = work_root / "manta" / sample_t
    strelka_dir = work_root / "strelka2" / sample_t
    merge_dir = strelka_dir / "merged"

    # dirs
    tool.judge_then_exec(run_sample_id, f"mkdir -p {manta_dir}", str(manta_dir))
    tool.judge_then_exec(run_sample_id, f"mkdir -p {strelka_dir}", str(strelka_dir))
    tool.judge_then_exec(run_sample_id, f"mkdir -p {merge_dir}", str(merge_dir))

    manta_runwf = str(manta_dir / "runWorkflow.py")
    strelka_runwf = str(strelka_dir / "runWorkflow.py")

    manta_cand = str(manta_dir / "results/variants/candidateSmallIndels.vcf.gz")
    snv_vcf = str(strelka_dir / "results/variants/somatic.snvs.vcf.gz")
    indel_vcf = str(strelka_dir / "results/variants/somatic.indels.vcf.gz")

    merged_vcf = str(merge_dir / f"{sample_t}.strelka.somatic.merged.vcf.gz")
    pass_vcf = str(merge_dir / f"{sample_t}.strelka.somatic.merged.PASS.vcf.gz")

    # 1) Manta configure
    tool.write_log("Strelka2: Manta configure", "info")
    cmd_manta_conf = (
        f"{args.config_manta} "
        f"--tumorBam {args.tumor_bam} "
        f"--normalBam {args.normal_bam} "
        f"--referenceFasta {args.ref} "
        f"--exome "
        f"--runDir {manta_dir}"
    )
    tool.judge_then_exec(run_sample_id, cmd_manta_conf, manta_runwf)

    # 2) Manta run
    tool.write_log("Strelka2: Manta run", "info")
    cmd_manta_run = f"{manta_runwf} -m local -j {args.threads}"
    tool.judge_then_exec(run_sample_id, cmd_manta_run, manta_cand)

    # 3) Strelka configure
    tool.write_log("Strelka2: configure", "info")
    cmd_strelka_conf = (
        f"{args.strelka_config} "
        f"--tumorBam {args.tumor_bam} "
        f"--normalBam {args.normal_bam} "
        f"--referenceFasta {args.ref} "
        f"--callRegions {args.call_regions} "
        f"--indelCandidates {manta_cand} "
        f"--exome "
        f"--runDir {strelka_dir}"
    )
    tool.judge_then_exec(run_sample_id, cmd_strelka_conf, strelka_runwf)

    # 4) Strelka run
    tool.write_log("Strelka2: run", "info")
    cmd_strelka_run = f"{strelka_runwf} -m local -j {args.threads}"
    tool.judge_then_exec(run_sample_id, cmd_strelka_run, snv_vcf)

    # 5) Merge
    tool.write_log("Strelka2: merge SNV/INDEL", "info")
    cmd_merge = (
        f"{args.bcftools} concat -a -Oz -o {merged_vcf} {snv_vcf} {indel_vcf} "
        f"&& {args.bcftools} index -t {merged_vcf}"
    )
    tool.judge_then_exec(run_sample_id, cmd_merge, merged_vcf)

    # 6) PASS
    tool.write_log("Strelka2: filter PASS", "info")
    cmd_pass = (
        f"{args.bcftools} view -f PASS {merged_vcf} -Oz -o {pass_vcf} "
        f"&& {args.bcftools} index -t {pass_vcf}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pass, pass_vcf)

    tool.write_log(f"Strelka2 done. PASS={pass_vcf}", "info")
    return {"pass_vcf": pass_vcf, "merged_vcf": merged_vcf, "snv_vcf": snv_vcf, "indel_vcf": indel_vcf}