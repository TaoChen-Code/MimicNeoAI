#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Example:
  python 03-salmon_quant.py \
    -s SAMPLE_T \
    -o ./salmon_quant \
    --ref-tx-fa novel/04a.trace_to_ref/ref.transcripts.fa \
    --contigs-fa novel/04a.trace_to_ref/contigs.annot.fa \
    --genome-fa /path/to/ref/GRCh38.p3.genome.fa \
    --fq1 /path/to/clean/SAMPLE_T.R1.QC.fq.gz \
    --fq2 /path/to/clean/SAMPLE_T.R2.QC.fq.gz \
    --threads 20 --kmer 31
"""

import argparse
import os
import sys
import hashlib
import subprocess
from pathlib import Path


def run(cmd):
    print("[CMD]", " ".join(map(str, cmd)), flush=True)
    subprocess.check_call(list(map(str, cmd)))


def check_file(p):
    if not Path(p).exists() or Path(p).stat().st_size == 0:
        sys.stderr.write(f"[ERR] missing or empty: {p}\n");
        sys.exit(1)


def _split_header(h):
    """Return (id_token, remaining_desc). Salmon uses the first whitespace-free token as target ID."""
    parts = h.split()
    if not parts:
        return h, ""
    return parts[0], " ".join(parts[1:])


def fasta_iter(fp):
    """Yield (raw_header_without_>, seq_string_no_newlines)."""
    h, seq = None, []
    with open(fp, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if h is not None:
            yield h, "".join(seq)


def write_merged_tx(ref_tx_fa, contigs_fa, out_fa,
                    dedup=True,
                    map_tsv=None,
                    annotate_dup=True,
                    max_dup_in_header=50):
    """
    Merge reference transcripts and contigs; optionally de-duplicate by sequence (md5).
    - If a contig has exactly the same sequence as a reference transcript, keep the ref entry,
      record the contig header as a duplicate in the kept record, and optionally write a mapping
      table (kept_id, dup_id, source) to map_tsv.
    - Salmon only takes the first token in the FASTA header as target_id; annotations are appended
      after a space and do not affect quant.sf target IDs.
    """
    seq_map = {}  # md5 -> kept_record

    # kept_record: {"source": "ref"/"contigs", "header": raw_header, "seq": "ATGC", "id": id_token, "dups": [(source,id,raw_header), ...]}

    def _ingest(fp, source):
        for raw_h, s in fasta_iter(fp):
            md5 = hashlib.md5(s.encode("utf-8")).hexdigest()
            tid, rest = _split_header(raw_h)
            rec = {"source": source, "header": raw_h, "seq": s, "id": tid, "dups": []}
            if not dedup:
                # No deduplication: treat every record independently
                seq_map[(md5, tid, source)] = rec
            else:
                if md5 not in seq_map:
                    seq_map[md5] = rec
                else:
                    # Sequence already seen: register as duplicate of the kept record
                    kept = seq_map[md5]
                    kept["dups"].append((source, tid, raw_h))

    # Ingest reference first, then contigs (so identical sequences prefer the ref record)
    _ingest(ref_tx_fa, "ref")
    _ingest(contigs_fa, "contigs")

    # Optional mapping table
    if map_tsv:
        with open(map_tsv, "w") as fo:
            fo.write("kept_id\tkept_source\tmd5\tdup_source\tdup_id\tdup_header\n")
            for key, kept in seq_map.items():
                # For dedup=False, key is (md5, tid, source); otherwise it's md5
                md5 = key if isinstance(key, str) else (key[0] if isinstance(key, tuple) else "NA")
                if isinstance(key, tuple) and len(key) == 3:
                    md5 = key[0]
                elif isinstance(key, str):
                    md5 = key
                else:
                    md5 = "NA"
                for (src, dup_id, dup_h) in kept.get("dups", []):
                    fo.write(f"{kept['id']}\t{kept['source']}\t{md5}\t{src}\t{dup_id}\t{dup_h}\n")
        print(f"[INFO] dedup map written: {map_tsv}", flush=True)

    # Write merged FASTA
    with open(out_fa, "w") as fo:
        items = seq_map.items()
        if dedup:
            # md5 -> kept_record
            items = [(k, v) for k, v in seq_map.items() if isinstance(k, str)]

        n_keep = 0
        for key, kept in items:
            tid, rest = _split_header(kept["header"])
            desc = rest
            if dedup and annotate_dup and kept["dups"]:
                # Append only contig-origin duplicates to avoid bloated headers from ref-internal dups
                contig_dup_ids = [dup_id for (src, dup_id, _raw) in kept["dups"] if src == "contigs"]
                if contig_dup_ids:
                    shown = contig_dup_ids[:max_dup_in_header]
                    overflow = len(contig_dup_ids) - len(shown)
                    tag = f"DUP_FROM_CONTIGS={','.join(shown)}"
                    if overflow > 0:
                        tag += f";DUP_TOTAL={len(contig_dup_ids)}"
                    desc = (desc + " " + tag).strip()

            # Salmon reads only the first token as ID; the trailing desc is annotation only
            fo.write(f">{tid}")
            if desc:
                fo.write(" " + desc)
            fo.write("\n")
            s = kept["seq"]
            for i in range(0, len(s), 60):
                fo.write(s[i:i + 60] + "\n")
            n_keep += 1

    print(f"[INFO] merged_tx written: kept={n_keep}, out={out_fa}", flush=True)


def write_decoys(genome_fa, decoys_txt):
    """Extract sequence names from genome FASTA headers (first token) into decoys.txt."""
    n = 0
    with open(genome_fa, "r") as fi, open(decoys_txt, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                tok = line[1:].strip().split()[0]
                fo.write(tok + "\n")
                n += 1
    print(f"[INFO] decoys.txt: {n} entries -> {decoys_txt}", flush=True)


def write_gentrome(merged_tx_fa, genome_fa, gentrome_fa):
    """gentrome = transcripts + genome concatenated in order."""
    with open(gentrome_fa, "w") as fo:
        for src in [merged_tx_fa, genome_fa]:
            with open(src, "r") as fi:
                for line in fi:
                    fo.write(line)
    print(f"[INFO] gentrome.fa written: {gentrome_fa}", flush=True)


def parse_args():
    ap = argparse.ArgumentParser(description="Build Salmon gentrome index (genome as decoys) and run quant.")
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("-o", "--outdir", required=True)

    # Inputs: reference transcripts + contigs + genome
    ap.add_argument("--ref-tx-fa", required=True, help="Reference transcript FASTA (coding + non-coding)")
    ap.add_argument("--contigs-fa", required=True, help="Trinity contigs from novel stage (or custom contigs)")
    ap.add_argument("--genome-fa", required=True, help="Reference genome FASTA (as decoys)")

    # Paired-end reads
    ap.add_argument("--fq1", required=True)
    ap.add_argument("--fq2", required=True)

    # Salmon options
    ap.add_argument("--threads", default="20")
    ap.add_argument("--kmer", default="31")
    ap.add_argument("--no-dedup", action="store_true", help="Disable sequence-level deduplication (enabled by default)")
    ap.add_argument("--no-annotate-dup", action="store_true", help="Do not append duplicate IDs in FASTA headers")
    ap.add_argument("--max-dup-in-header", type=int, default=50,
                    help="Max number of duplicate IDs to append in header (default: 50)")

    return ap.parse_args()


def main():
    a = parse_args()
    outdir = Path(a.outdir);
    outdir.mkdir(parents=True, exist_ok=True)

    # Basic checks
    for p in [a.ref_tx_fa, a.contigs_fa, a.genome_fa, a.fq1, a.fq2]:
        check_file(p)

    merged_tx = outdir / f"{a.sample}.merged_tx.fa"
    dedup_map = outdir / f"{a.sample}.merged_tx.dedup.map.tsv"
    gentrome = outdir / f"{a.sample}.gentrome.fa"
    decoys = outdir / f"{a.sample}.decoys.txt"
    idx_dir = outdir / "salmon_index"
    quant_dir = outdir / "salmon_quant"

    # 1) Merge transcripts (with optional deduplication)
    if not merged_tx.exists():
        write_merged_tx(
            a.ref_tx_fa, a.contigs_fa, merged_tx,
            dedup=(not a.no_dedup),
            map_tsv=str(dedup_map),
            annotate_dup=(not a.no_annotate_dup),
            max_dup_in_header=a.max_dup_in_header
        )
    else:
        print(f"[SKIP] merged_tx exists: {merged_tx}")

    # 2) decoys.txt from genome
    if not decoys.exists():
        write_decoys(a.genome_fa, decoys)
    else:
        print(f"[SKIP] decoys exists: {decoys}")

    # 3) gentrome.fa = merged_tx + genome
    if not gentrome.exists():
        write_gentrome(merged_tx, a.genome_fa, gentrome)
    else:
        print(f"[SKIP] gentrome exists: {gentrome}")

    # 4) salmon index (gentrome + decoys)
    if not idx_dir.exists():
        run(["salmon", "index",
             "-t", str(gentrome),
             "-d", str(decoys),
             "-i", str(idx_dir),
             "-k", str(a.kmer),
             "-p", str(a.threads)])
    else:
        print(f"[SKIP] index exists: {idx_dir}")

    # 5) salmon quant (paired-end reads, library type auto-detection + validateMappings)
    quant_dir.mkdir(exist_ok=True)
    quant_sf = quant_dir / "quant.sf"

    if quant_sf.exists() and quant_sf.stat().st_size > 0:
        print(f"[SKIP] quant exists: {quant_sf}")
    else:
        run([
            "salmon", "quant",
            "-i", str(idx_dir),
            "-l", "A",
            "-1", a.fq1, "-2", a.fq2,
            "-p", str(a.threads),
            "--validateMappings",
            "--gcBias",
            "--seqBias",
            "-o", str(quant_dir)
        ])

    print("[DONE] quant:", quant_dir)


if __name__ == "__main__":
    main()
