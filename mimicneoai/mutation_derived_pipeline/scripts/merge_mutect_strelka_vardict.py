#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
merge_mutect_strelka_vardict.py

Purpose
-------
Normalize 3 PASS VCFs (Mutect2 / Strelka2 / VarDict) and generate:
  01.norm/              sorted+normalized VCFs (bgz + tabix)
  02.shared_mutect/     intersection VCF (>=2 callers; records written from Mutect)
  03.only/              caller-specific-only VCFs (by position)
  04.bed/               BED with tumor/normal AF/AD + pad200 BED + headers.txt

Design notes
------------
- This file is intended to be IMPORTED and called by pipeline code:
    from ...merge_mutect_strelka_vardict import merge_mutect_strelka_vardict
  so that ALL external commands are executed via tool.judge_then_exec(...).
- We do NOT rely on capturing stdout from tool.judge_then_exec. For bcftools query,
  we redirect stdout to a file and parse it in Python.
"""

import os
import pathlib
from typing import Optional, Dict, Any


def _must_exist(path: str, label: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"[ERROR] {label} not found: {path}")


def must_be_pass_vcfgz(path: str, label: str) -> None:
    """Guard: ensure the input looks like PASS-filtered bgzipped VCF."""
    _must_exist(path, label)
    base = os.path.basename(path)
    if not base.endswith(".vcf.gz"):
        raise ValueError(f"[ERROR] {label} must be .vcf.gz: {path}")
    if ("PASS" not in base) and (".pass" not in base.lower()):
        raise ValueError(f"[ERROR] {label} does not look like PASS-filtered VCF: {path}")
    # index is recommended but not mandatory here
    # (downstream steps tabix/index during normalization anyway)


def merge_mutect_strelka_vardict(
    *,
    tool,
    run_sample_id: str,
    ref: str,
    mutect: str,
    strelka: str,
    vardict: str,
    tumor: str,
    normal: str,
    out_dir: str,
    prefix: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Normalize 3 PASS VCFs (Mutect/Strelka/VarDict) -> intersections/only -> export AF/AD BED (+pad200).

    Parameters
    ----------
    tool : pipeline tool object with judge_then_exec/run logging capabilities
    run_sample_id : used by tool for logging/de-dup
    ref : reference FASTA
    mutect/strelka/vardict : PASS-filtered bgzipped VCFs (.vcf.gz)
    tumor/normal : sample names used in VCF FORMAT fields (bcftools query -s)
    out_dir : output root directory (will create subdirs)
    prefix : output prefix (default=tumor)

    Returns
    -------
    Dict with key output paths.
    """
    if tool is None:
        raise ValueError("[ERROR] tool is required")
    if not run_sample_id:
        raise ValueError("[ERROR] run_sample_id is required")

    _must_exist(ref, "ref")
    must_be_pass_vcfgz(mutect, "mutect")
    must_be_pass_vcfgz(strelka, "strelka")
    must_be_pass_vcfgz(vardict, "vardict")

    out_root = pathlib.Path(out_dir)
    tool.judge_then_exec(run_sample_id, f"mkdir -p {out_root}", str(out_root))

    prefix = prefix or tumor

    # subdirs
    d_norm = out_root / "01.norm"
    d_shared = out_root / "02.shared_mutect"
    d_only = out_root / "03.only"
    d_bed = out_root / "04.bed"
    for d in (d_norm, d_shared, d_only, d_bed):
        tool.judge_then_exec(run_sample_id, f"mkdir -p {d}", str(d))

    # basenames
    M_BASE = f"{prefix}.mutect"
    S_BASE = f"{prefix}.strelka"
    V_BASE = f"{prefix}.vardict"

    # ---- 01. normalize (sort -> norm split) ----
    def do_norm(in_vcfgz: str, base: str) -> str:
        sorted_vcf = d_norm / f"{base}.sorted.vcf.gz"
        norm_vcf = d_norm / f"{base}.norm.vcf.gz"

        tool.write_log(f"Merge3Callers: normalize {os.path.basename(in_vcfgz)}", "info")

        cmd_sort = f"bcftools sort -Oz -o {sorted_vcf} {in_vcfgz} && tabix -p vcf {sorted_vcf}"
        tool.judge_then_exec(run_sample_id, cmd_sort, str(sorted_vcf))

        cmd_norm = f"bcftools norm -f {ref} -m -both -Oz -o {norm_vcf} {sorted_vcf} && tabix -p vcf {norm_vcf}"
        tool.judge_then_exec(run_sample_id, cmd_norm, str(norm_vcf))

        return str(norm_vcf)

    MN = do_norm(mutect, M_BASE)
    SN = do_norm(strelka, S_BASE)
    VN = do_norm(vardict, V_BASE)

    # ---- 02. intersections & only (by position) ----
    shared_vcf = d_shared / f"{prefix}.shared.vcf.gz"
    only_m = d_only / f"{prefix}.mutect.only.vcf.gz"
    only_s = d_only / f"{prefix}.strelka.only.vcf.gz"
    only_v = d_only / f"{prefix}.vardict.only.vcf.gz"

    tool.write_log("Merge3Callers: bcftools isec (shared >=2 callers; write records from Mutect)", "info")
    cmd_shared = (
        f"bcftools isec -c none -n+2 -w1 {MN} {SN} {VN} -O z -o {shared_vcf} "
        f"&& tabix -p vcf {shared_vcf}"
    )
    tool.judge_then_exec(run_sample_id, cmd_shared, str(shared_vcf))

    tool.write_log("Merge3Callers: bcftools isec (only sets)", "info")
    cmd_only_m = f"bcftools isec -c none -n=1 -w1 {MN} {SN} {VN} -O z -o {only_m} && tabix -p vcf {only_m}"
    cmd_only_s = f"bcftools isec -c none -n=1 -w2 {MN} {SN} {VN} -O z -o {only_s} && tabix -p vcf {only_s}"
    cmd_only_v = f"bcftools isec -c none -n=1 -w3 {MN} {SN} {VN} -O z -o {only_v} && tabix -p vcf {only_v}"
    tool.judge_then_exec(run_sample_id, cmd_only_m, str(only_m))
    tool.judge_then_exec(run_sample_id, cmd_only_s, str(only_s))
    tool.judge_then_exec(run_sample_id, cmd_only_v, str(only_v))

    # ---- 03. export AF+AD BED from shared (Mutect-based) ----
    tool.write_log("Merge3Callers: export AF/AD BED", "info")

    query_tsv = d_bed / f"{prefix}.shared.AF_AD.tsv"
    cmd_query = (
        f"bcftools query -s {tumor},{normal} "
        f"-f '%CHROM\\t%POS\\t[%AF\\t]\\t[%AD\\t]\\n' "
        f"{shared_vcf} > {query_tsv}"
    )
    tool.judge_then_exec(run_sample_id, cmd_query, str(query_tsv))

    out_bed = d_bed / f"{prefix}.shared.AF_AD.bed"
    out_bed_pad = d_bed / f"{prefix}.shared.AF_AD.pad200.bed"
    headers_txt = d_bed / "headers.txt"

    def split_ad(x: str):
        if x and ("," in x):
            a = x.split(",", 1)
            return a[0], a[1]
        return ".", "."

    # TSV -> BED
    with open(query_tsv, "r") as fin, open(out_bed, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")

            # expected: CHROM POS T_AF N_AF T_AD N_AD
            chrom = parts[0]
            pos1 = int(parts[1])          # 1-based
            start0 = max(0, pos1 - 1)     # 0-based

            t_af = parts[2] if len(parts) > 2 and parts[2] else "."
            n_af = parts[3] if len(parts) > 3 and parts[3] else "."
            t_ad = parts[4] if len(parts) > 4 and parts[4] else "."
            n_ad = parts[5] if len(parts) > 5 and parts[5] else "."

            ta_ref, ta_alt = split_ad(t_ad)
            na_ref, na_alt = split_ad(n_ad)

            fout.write(
                f"{chrom}\t{start0}\t{pos1}\t{t_af}\t{n_af}\t{ta_ref}\t{ta_alt}\t{na_ref}\t{na_alt}\n"
            )

    # pad200
    with open(out_bed, "r") as fin, open(out_bed_pad, "w") as fout:
        for l in fin:
            chrom, start0, end_pos, t_af, n_af, ta_ref, ta_alt, na_ref, na_alt = l.rstrip("\n").split("\t")
            start0 = int(start0)
            end_pos = int(end_pos)
            pad_start = max(0, start0 - 200)
            pad_end = end_pos + 200
            fout.write(
                f"{chrom}\t{pad_start}\t{pad_end}\t{end_pos}\t{t_af}\t{n_af}\t{ta_ref}\t{ta_alt}\t{na_ref}\t{na_alt}\n"
            )

    # headers
    with open(headers_txt, "w") as fh:
        fh.write(
            "shared.AF_AD.bed\tchrom\tstart_0based\tend_pos(=POS)\t"
            "tumor_AF\tnormal_AF\ttumor_AD_ref\ttumor_AD_alt\tnormal_AD_ref\tnormal_AD_alt\n"
        )
        fh.write(
            "shared.AF_AD.pad200.bed\tchrom\tpad_start_0based\tpad_end_exclusive\torig_pos_1based\t"
            "tumor_AF\tnormal_AF\ttumor_AD_ref\ttumor_AD_alt\tnormal_AD_ref\tnormal_AD_alt\n"
        )

    tool.write_log(f"Merge3Callers done. Outputs in: {out_root}", "info")
    return {
        "out_dir": str(out_root),
        "norm_mutect": MN,
        "norm_strelka": SN,
        "norm_vardict": VN,
        "shared_vcf": str(shared_vcf),
        "only_mutect_vcf": str(only_m),
        "only_strelka_vcf": str(only_s),
        "only_vardict_vcf": str(only_v),
        "query_tsv": str(query_tsv),
        "bed": str(out_bed),
        "bed_pad200": str(out_bed_pad),
        "headers": str(headers_txt),
    }


# ------------------------- Optional CLI (not used by pipeline) -------------------------
def main(argv=None):
    """
    Optional CLI wrapper.
    NOTE: Your pipeline should call merge_mutect_strelka_vardict(...) directly with tool.
    """
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        description="(CLI wrapper) This script is designed to be used inside the pipeline with tool."
    )
    ap.add_argument("--ref", required=True)
    ap.add_argument("--mutect", required=True)
    ap.add_argument("--strelka", required=True)
    ap.add_argument("--vardict", required=True)
    ap.add_argument("--tumor", required=True)
    ap.add_argument("--normal", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default=None)
    _ = ap.parse_args(argv)

    sys.stderr.write(
        "[ERROR] This script is intended to be called from pipeline with a tool object.\n"
        "Import and call merge_mutect_strelka_vardict(...).\n"
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())