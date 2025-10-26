#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example:
python 04-salmon_quant_control.py \
  --fq1 /path/to/control_R1.fastq.gz \
  --fq2 /path/to/control_R2.fastq.gz \
  -i /path/to/salmon_index \
  -o /path/to/output_dir \
  -p 20
"""

import argparse
import sys
import subprocess
from pathlib import Path


def run(cmd):
    """Execute a shell command and print it for reproducibility."""
    print("[CMD]", " ".join(map(str, cmd)), flush=True)
    subprocess.check_call(list(map(str, cmd)))


def check_path_exists(p):
    """Exit with an error if a required path does not exist."""
    if not Path(p).exists():
        sys.stderr.write(f"[ERR] Not found: {p}\n")
        sys.exit(1)


def parse_args():
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Run salmon quant on control FASTQs using an existing index."
    )
    ap.add_argument("--fq1", required=True, help="Control R1 FASTQ (gz supported)")
    ap.add_argument("--fq2", required=True, help="Control R2 FASTQ (gz supported)")
    ap.add_argument(
        "-i", "--index",
        default="salmon_index",
        help="Path to an existing salmon index (default: salmon_index)"
    )
    ap.add_argument(
        "-o", "--outdir",
        default="salmon_quant/salmon_quant_control",
        help="Output directory (default: salmon_quant/salmon_quant_control)"
    )
    ap.add_argument(
        "-p", "--threads",
        default="20",
        help="Number of threads (default: 20)"
    )
    return ap.parse_args()


def main():
    a = parse_args()

    # Basic existence checks
    check_path_exists(a.index)
    check_path_exists(a.fq1)
    check_path_exists(a.fq2)

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    quant_sf = outdir / "quant.sf"

    # Skip if a previous successful run is detected
    if quant_sf.exists() and quant_sf.stat().st_size > 0:
        print(f"[SKIP] Detected existing result: {quant_sf}")
        print("[INFO] To re-run, remove existing outputs or specify a new --outdir.")
        return

    # Quantification with automatic library type and common bias corrections
    run([
        "salmon", "quant",
        "-i", str(a.index),
        "-l", "A",
        "-1", a.fq1, "-2", a.fq2,
        "-p", str(a.threads),
        "--validateMappings",
        "--gcBias",
        "--seqBias",
        "-o", str(outdir)
    ])

    print("[DONE] Quantification completed ->", outdir)
    print("  - quant.sf")
    print("  - aux_info/meta_info.json")
    print("  - logs/salmon_quant.log")


if __name__ == "__main__":
    main()
