#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example:
  python 01-alignment.py \
    -s Sample_001_T \
    --clean-dir ./00-clean \
    --genome-dir star_and_rsem/GRCh38_gencode.v23 \
    --out-root ./01-star \
    -p 20
"""

import argparse
import sys
import subprocess
from pathlib import Path

def run(cmd, shell=False):
    """Run a command and print it; raise on non-zero exit."""
    print("[CMD]", cmd if shell else " ".join(map(str, cmd)), flush=True)
    if shell:
        subprocess.check_call(cmd, shell=True, executable="/bin/bash")
    else:
        subprocess.check_call(cmd)

def parse_args():
    ap = argparse.ArgumentParser(
        description="STAR alignment using cleaned FASTQs produced by 00-QC.py"
    )
    ap.add_argument("-s", "--sample", required=True, help="Sample name, e.g., SAMPLE-T-RNA")
    ap.add_argument("--clean-dir", default="./00-clean", help="Directory of cleaned FASTQs (*.QC.fq.gz)")
    # ap.add_argument("--star-bin", default="STAR", help="Path to STAR binary (default: STAR in PATH)")
    ap.add_argument("--genome-dir", required=True, help="STAR --genomeDir (RSEM/STAR index)")
    ap.add_argument("-p", "--threads", default="20", help="Number of threads for STAR")
    # Output root; keep layout consistent with downstream scripts: 01.star/<sample>/<sample>.star/
    ap.add_argument("--out-root", default="./01.star", help="Output root directory (default: ./01.star)")
    return ap.parse_args()

def main():
    a = parse_args()

    clean_dir = Path(a.clean_dir).resolve()
    r1 = clean_dir / f"{a.sample}.R1.QC.fq.gz"
    r2 = clean_dir / f"{a.sample}.R2.QC.fq.gz"
    if not (r1.exists() and r2.exists()):
        sys.stderr.write(f"[ERR] Cleaned FASTQs not found:\n  {r1}\n  {r2}\n")
        sys.stderr.write("      (Run 00-QC.py first, or provide correct --clean-dir / --sample)\n")
        sys.exit(1)

    # STAR output directory and filename prefix (aligned with expectations in downstream pipeline)
    star_dir = Path(a.out_root, a.sample, f"{a.sample}.star").resolve()
    star_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = str(star_dir / a.sample)  # STAR will produce <prefix>Aligned.out.bam
    aligned_bam = star_dir / f"{a.sample}Aligned.out.bam"

    # Skip if the expected BAM already exists and is non-empty
    if aligned_bam.exists() and aligned_bam.stat().st_size > 0:
        print(f"[SKIP] Found existing BAM: {aligned_bam}")
    else:
        cmd_star = fr"""STAR \
         --genomeDir {a.genome_dir} \
         --outSAMunmapped Within \
         --outFilterType BySJout \
         --outSAMattributes NH HI AS NM MD \
         --outFilterMultimapNmax 20 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --sjdbScore 1 \
         --runThreadN {a.threads} \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 40000000000 \
         --outSAMtype BAM Unsorted \
         --quantMode TranscriptomeSAM \
         --outSAMheaderHD \@HD VN:1.4 SO:unsorted \
         --outFileNamePrefix {out_prefix} \
         --readFilesCommand zcat \
         --readFilesIn {r1} {r2}
        """
        run(cmd_star, shell=True)

    if aligned_bam.exists() and aligned_bam.stat().st_size > 0:
        print("[DONE] STAR finished.")
        print(" BAM   :", aligned_bam)
        print(" Prefix outputs under:", star_dir)
    else:
        sys.stderr.write("[ERR] STAR did not produce Aligned.out.bam as expected.\n")
        sys.exit(2)

if __name__ == "__main__":
    main()
