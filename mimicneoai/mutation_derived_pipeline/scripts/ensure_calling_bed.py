import os
import shutil
import subprocess
from typing import Tuple


def _run(cmd: list, tool=None, run_sample_id: str = "", check: bool = True) -> None:
    """Run a command, raise on failure."""
    if tool is not None:
        tool.write_log(f"CMD: {' '.join(cmd)}", "info")
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and p.returncode != 0:
        msg = f"Command failed ({p.returncode}): {' '.join(cmd)}\nSTDERR:\n{p.stderr}"
        if tool is not None:
            tool.write_log(msg, "error")
        raise RuntimeError(msg)


def _ensure_calling_bed(
    run_sample_id: str,
    tool,
    dir_varcall_root_tumor: str,
    bed_file: str,
    pad_bp: int,
    bcftools: str = "bcftools",
) -> Tuple[str, str]:
    """
    Create padded calling BED files under:
      <dir_varcall_root_tumor>/bed/

    Outputs:
      calling_bed:     .../bed/calling.bed
      call_regions_gz: .../bed/calling.bed.gz   (bgzip)
      index:           .../bed/calling.bed.gz.tbi  (tabix -p bed)

    Returns:
      (calling_bed, call_regions_bed_gz)
    """
    bed_file = str(bed_file).strip()
    if not bed_file:
        raise ValueError("bed_file is empty")
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"bed_file not found: {bed_file}")

    pad_bp = int(pad_bp)
    if pad_bp < 0:
        raise ValueError(f"pad_bp must be >=0, got {pad_bp}")

    bed_dir = os.path.join(dir_varcall_root_tumor, "bed")
    os.makedirs(bed_dir, exist_ok=True)

    calling_bed = os.path.join(bed_dir, f"calling.pad{pad_bp}.bed")
    calling_bed_gz = calling_bed + ".gz"
    calling_bed_tbi = calling_bed_gz + ".tbi"

    # If already prepared, reuse
    if os.path.exists(calling_bed) and os.path.exists(calling_bed_gz) and os.path.exists(calling_bed_tbi):
        tool.write_log(f"Reuse calling BED: {calling_bed_gz}", "info")
        return calling_bed, calling_bed_gz

    tool.write_log(f"Prepare calling BED with pad_bp={pad_bp}: {bed_file}", "info")

    # ---- 1) Read + pad ----
    tmp_unsorted = os.path.join(bed_dir, f".tmp.calling.pad{pad_bp}.unsorted.bed")
    with open(bed_file, "r") as fin, open(tmp_unsorted, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except Exception:
                continue

            # pad
            start2 = max(0, start - pad_bp)
            end2 = end + pad_bp

            # keep remaining columns if any
            rest = cols[3:]
            out_cols = [chrom, str(start2), str(end2)] + rest
            fout.write("\t".join(out_cols) + "\n")

    # ---- 2) Sort (required for tabix) ----
    # Prefer system sort: sort -k1,1 -k2,2n
    # (If you need "chr" aware ordering, we can swap to bedtools sort with genome file.)
    tmp_sorted = os.path.join(bed_dir, f".tmp.calling.pad{pad_bp}.sorted.bed")
    _run(["bash", "-lc", f"sort -k1,1 -k2,2n {tmp_unsorted} > {tmp_sorted}"], tool, run_sample_id)

    # atomic move
    os.replace(tmp_sorted, calling_bed)
    try:
        os.remove(tmp_unsorted)
    except OSError:
        pass

    # ---- 3) bgzip + tabix ----
    # We need bgzip+tabix (preferred), otherwise try bcftools (only if it provides bgzip/tabix wrappers)
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")

    if bgzip and tabix:
        _run(["bash", "-lc", f"bgzip -c {calling_bed} > {calling_bed_gz}"], tool, run_sample_id)
        _run(["tabix", "-f", "-p", "bed", calling_bed_gz], tool, run_sample_id)
    else:
        # fallback: if bcftools exists and provides bgzip/tabix (some installs do not)
        # 1) try `bcftools view` won't help; need bgzip specifically.
        # We try `bcftools`-bundled `bgzip/tabix` if accessible via PATH; otherwise fail clearly.
        tool.write_log(
            "bgzip/tabix not found in PATH. Please install htslib (bgzip/tabix) or ensure they are in PATH.",
            "error",
        )
        raise RuntimeError("Missing bgzip/tabix; cannot create call_regions_bed_gz with .tbi index")

    if not (os.path.exists(calling_bed_gz) and os.path.exists(calling_bed_tbi)):
        raise RuntimeError(f"Failed to generate {calling_bed_gz} and/or {calling_bed_tbi}")

    tool.write_log(f"Calling BED ready: {calling_bed_gz}", "info")
    return calling_bed, calling_bed_gz