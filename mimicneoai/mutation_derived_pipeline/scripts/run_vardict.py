import argparse
import os
import pathlib
from typing import Optional


def vardict_main(argv=None, tool=None, run_sample_id: Optional[str] = None):
    """
    VarDict paired:
      vardict-java -> testsomatic.R -> var2vcf_paired.pl -> bcftools PASS -> index
    All commands executed via tool.judge_then_exec (tool handles logs).
    """
    if tool is None:
        raise ValueError("vardict_main requires tool")

    ap = argparse.ArgumentParser(
        description="VarDict paired + testsomatic.R + var2vcf_paired.pl + bcftools PASS (tool version)"
    )
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--bed", required=True, help="Calling BED (.bed or .bed.gz)")
    ap.add_argument("--tumor-bam", required=True, help="Tumor BAM (BQSR)")
    ap.add_argument("--normal-bam", required=True, help="Normal BAM (BQSR)")
    ap.add_argument("--sample-t", required=True, help="Tumor sample name")
    ap.add_argument("--sample-n", required=True, help="Normal sample name")
    ap.add_argument("--af-thr", default="0.02", help="Min AF (default 0.02)")
    ap.add_argument("--threads", type=int, default=10, help="Threads (default 10)")
    ap.add_argument("--out-dir", default=".", help="Output dir")
    ap.add_argument("--prefix", default=None, help="Output prefix (default sample-t)")
    ap.add_argument("--bcftools", default="bcftools", help="bcftools path/name")
    args = ap.parse_args(argv)

    if not run_sample_id:
        run_sample_id = f"{args.sample_t},{args.sample_n}"

    out_dir = pathlib.Path(args.out_dir)
    tool.judge_then_exec(run_sample_id, f"mkdir -p {out_dir}", str(out_dir))

    prefix = args.prefix if args.prefix else args.sample_t
    out_vcf = str(out_dir / f"{prefix}.vardict.paired.vcf")
    out_pass_vcfgz = str(out_dir / f"{prefix}.vardict.paired.PASS.vcf.gz")

    # 1) VarDict pipeline -> VCF
    tool.write_log("VarDict: vardict-java | testsomatic.R | var2vcf_paired.pl", "info")
    cmd_pipeline = (
        f"vardict-java "
        f"-th {args.threads} "
        f"-G {args.ref} "
        f"-f {args.af_thr} "
        f"-N {args.sample_t} "
        f"-b '{args.tumor_bam}|{args.normal_bam}' "
        f"-z -c 1 -S 2 -E 3 -g 4 {args.bed} "
        f"| testsomatic.R "
        f"| var2vcf_paired.pl -N '{args.sample_t}|{args.sample_n}' -f {args.af_thr} "
        f"> {out_vcf}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pipeline, out_vcf)

    # 2) PASS only + index
    tool.write_log("VarDict: bcftools PASS + index", "info")
    cmd_pass = (
        f"{args.bcftools} view -f PASS -Oz -o {out_pass_vcfgz} {out_vcf} "
        f"&& {args.bcftools} index -t {out_pass_vcfgz}"
    )
    tool.judge_then_exec(run_sample_id, cmd_pass, out_pass_vcfgz)

    tool.write_log(f"VarDict done. PASS={out_pass_vcfgz}", "info")
    return {"vcf": out_vcf, "pass_vcfgz": out_pass_vcfgz}