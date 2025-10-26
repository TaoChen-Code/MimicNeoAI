#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example (sanitized):
python 05-hla_typing.py \
  -s SAMPLE_ID \
  --r1 /path/to/clean_reads/SAMPLE_ID.R1.QC.fq.gz \
  --r2 /path/to/clean_reads/SAMPLE_ID.R2.QC.fq.gz \
  --output-dir /path/to/project/results \
  --step-name-hla 05-hla_typing_test \
  -t 20 \
  --freq-data-dir /path/to/hlahd/freq_data/ \
  --HLA-gene      /path/to/hlahd/HLA_gene.split.txt \
  --dictionary    /path/to/hlahd/dictionary/ \
  --hla-gen       /path/to/hla_gen/hla_gen
"""

import os
import argparse

def run_cmd(cmd: str):
    """Run a shell command; raise if exit code != 0."""
    print("[CMD]", cmd, flush=True)
    rc = os.system(cmd)
    if rc != 0:
        raise RuntimeError(f"Command failed (exit {rc}): {cmd}")

def ensure_dir(p: str):
    """Create directory if it doesn't exist."""
    os.makedirs(p, exist_ok=True)

def is_directory_empty(dir_path: str) -> bool:
    """Return True if directory does not exist or has no files."""
    return (not os.path.isdir(dir_path)) or (not os.listdir(dir_path))

def build_argparser():
    ap = argparse.ArgumentParser(
        description="HLA typing pipeline using bowtie2/samtools/awk/HLA-HD (invoked via os.system)."
    )
    ap.add_argument("-s", "--sample", required=True, help="Sample ID")
    ap.add_argument("--r1", required=True, help="R1 FASTQ(.gz)")
    ap.add_argument("--r2", required=True, help="R2 FASTQ(.gz)")
    ap.add_argument("--output-dir", required=True, help="Output root directory")
    ap.add_argument("--step-name-hla", default="05-hla_typing", help="Step folder name (default: 05-hla_typing)")
    ap.add_argument("-t", "--thread", type=int, default=16, help="Number of threads")

    # HLA-HD resources
    ap.add_argument("--freq-data-dir", required=True, help="HLA-HD frequency data directory")
    ap.add_argument("--HLA-gene",      required=True, help="HLA-HD HLA_gene file")
    ap.add_argument("--dictionary",    required=True, help="HLA-HD dictionary directory")
    ap.add_argument("--hla-gen",       required=True, help="Bowtie2 index prefix (-x)")
    return ap

def main():
    args = build_argparser().parse_args()
    sample      = args.sample
    r1          = args.r1
    r2          = args.r2
    output_path = args.output_dir.rstrip("/") + "/"
    step_name   = args.step_name_hla
    thread      = int(args.thread)

    freq_data_dir = args.freq_data_dir
    HLA_gene      = args.HLA_gene
    dictionary    = args.dictionary
    hla_gen       = args.hla_gen

    # Output layout
    output_hla = f"{output_path}/{step_name}/"
    fastq_dir  = f"{output_hla}/{sample}/fastq/"
    ensure_dir(fastq_dir)

    # 1) Map to HLA reference with bowtie2 → keep mapped pairs with samtools → convert back to FASTQ
    #    Note: some bowtie2 versions use -S (stdout SAM) instead of -o (output path). Detect and fallback.
    cmd_1 = f"bowtie2 -p {thread} -x {hla_gen} -1 {r1} -2 {r2} -o {fastq_dir}/{sample}.hlamap.sam"
    if os.system("bowtie2 -h 2>&1 | grep -q ' -o '") != 0:
        cmd_1 = f"bowtie2 -p {thread} -x {hla_gen} -1 {r1} -2 {r2} -S {fastq_dir}/{sample}.hlamap.sam"

    cmd_2 = f"samtools view -@ {thread} -h -F 4 {fastq_dir}/{sample}.hlamap.sam > {fastq_dir}/{sample}.mapped.sam"
    cmd_3 = (
        f"samtools fastq -@ {thread} "
        f"-1 {fastq_dir}/{sample}.hlatmp.1.fastq "
        f"-2 {fastq_dir}/{sample}.hlatmp.2.fastq "
        f"-0 /dev/null -s /dev/null -n {fastq_dir}/{sample}.mapped.sam"
    )

    # 2) Normalize read headers for HLA-HD (replace '/1' → ' 1', '/2' → ' 2' in FASTQ headers)
    # awk_cmd_1 = r"'{if(NR%4 == 1){O=$0;gsub(\"/1\",\" 1\",O);print O}else{print $0}}'"
    # awk_cmd_2 = r"'{if(NR%4 == 1){O=$0;gsub(\"/2\",\" 2\",O);print O}else{print $0}}'"
    # cmd_4 = f"cat {fastq_dir}/{sample}.hlatmp.1.fastq | awk {awk_cmd_1} > {fastq_dir}/{sample}.hla.1.fastq"
    # cmd_5 = f"cat {fastq_dir}/{sample}.hlatmp.2.fastq | awk {awk_cmd_2} > {fastq_dir}/{sample}.hla.2.fastq"

    # FASTQ processing commands
    awk_cmd_1 = "'{if(NR%4 == 1){O=$0;gsub(\"/1\",\" 1\",O);print O}else{print $0}}'"
    awk_cmd_2 = "'{if(NR%4 == 1){O=$0;gsub(\"/2\",\" 2\",O);print O}else{print $0}}'"
    cmd_4 = f"cat {fastq_dir}/{sample}.hlatmp.1.fastq |awk {awk_cmd_1} > {fastq_dir}/{sample}.hla.1.fastq"
    cmd_5 = f"cat {fastq_dir}/{sample}.hlatmp.2.fastq |awk {awk_cmd_2} > {fastq_dir}/{sample}.hla.2.fastq"

    # 3) Run HLA-HD main workflow
    cmd_6 = (
        f"hlahd.sh -t {thread} -m 30 -f {freq_data_dir} "
        f"{fastq_dir}/{sample}.hla.1.fastq {fastq_dir}/{sample}.hla.2.fastq "
        f"{HLA_gene} {dictionary} {sample} {output_hla}"
    )
    final_result = f"{output_hla}/{sample}/result/{sample}_final.result.txt"

    # Execute only if final result doesn't exist
    if not os.path.exists(final_result):
        run_cmd(cmd_1)
        run_cmd(cmd_2)
        run_cmd(cmd_3)
        run_cmd(cmd_4)
        run_cmd(cmd_5)
        run_cmd(cmd_6)
    else:
        print(f"[INFO] Skip HLA-HD (exists): {final_result}")

    # 4) Cleanup temporary files if present
    mapfile_dir = f"{output_hla}/{sample}/mapfile/"
    intron_dir  = f"{output_hla}/{sample}/intron/"
    exon_dir    = f"{output_hla}/{sample}/exon/"

    if os.path.exists(f'{fastq_dir}/{sample}.hlamap.sam'):
        run_cmd(f"rm -f {fastq_dir}/*")
    if os.path.exists(mapfile_dir) and (not is_directory_empty(mapfile_dir)):
        run_cmd(f"rm -f {mapfile_dir}/*")
    if os.path.exists(intron_dir) and (not is_directory_empty(intron_dir)):
        run_cmd(f"rm -f {intron_dir}/*")
    if os.path.exists(exon_dir) and (not is_directory_empty(exon_dir)):
        run_cmd(f"rm -f {exon_dir}/*")

    print("[DONE] HLA typing finished.", flush=True)

if __name__ == "__main__":
    main()
