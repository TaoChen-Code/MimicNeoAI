#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gzip
import argparse
from pathlib import Path
from typing import Optional, Set
import pandas as pd
import numpy as np

"""
Example:
python 06-sORF-encoded_peptides.py \
  --quant-tumor /path/to/tumor/quant.sf \
  --quant-ctrl  /path/to/control/quant.sf \
  --dedup-map   /path/to/merged_tx.dedup.map.tsv \
  --known-tx-fa /path/to/lncRNA.consensus.fa \
  --pep-known   /path/to/known.lncRNA.sORFs.pep \
  --pep-novel   /path/to/novel.lncRNA.sORFs.pep \
  --out-dir     /path/to/output_dir \
  --out-name    sample.aeSEPs.pep \
  --save-table  sample.aberrant_noncoding.annot.csv \
  --min-tpm-tumor 1.0 \
  --max-tpm-ctrl 0.5 \
  --min-log2fc 2.0
"""


# -------------------------
# Small utilities
# -------------------------
def check_file(p: str):
    """Raise if a required file is missing or empty."""
    if not Path(p).exists() or Path(p).stat().st_size == 0:
        raise FileNotFoundError(f"[ERR] missing or empty: {p}")


def fasta_iter(fp):
    """Yield (header_without_>, sequence_no_newlines) for a (possibly gzipped) FASTA."""
    h, seq = None, []
    op = gzip.open if str(fp).endswith(".gz") else open
    with op(fp, "rt") as fh:
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


def _base_name_from_header(h: str) -> str:
    """Take the first token of a header and strip the trailing '.pX' suffix if present."""
    tok = h.split()[0]
    return tok.rsplit(".p", 1)[0]


def load_fasta_ids(fa_path: str) -> Set[str]:
    """
    Read a (possibly gzipped) FASTA and return a set of the first token in each header
    (e.g., ENST...).
    """
    ids: Set[str] = set()
    op = gzip.open if str(fa_path).endswith(".gz") else open
    with op(fa_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                tid = line[1:].strip().split()[0]
                if tid:
                    ids.add(tid)
    return ids


# -------------------------
# Core functions
# -------------------------
def annotate_with_dedup_map(res_df: pd.DataFrame, dedup_map_tsv: str) -> pd.DataFrame:
    """
    Merge de-duplication mapping info into res_df (expects a 'Name' column).

    Adds:
      - kept_source: source of the kept entry (ref/contigs)
      - has_trinity_dup, n_trinity_dup, trinity_list: collapsed duplicates from contigs/TRINITY
      - has_ref_dup,     n_ref_dup,     ref_dup_list: collapsed duplicates from reference (same seq other ENSTs)
    """
    m = pd.read_csv(dedup_map_tsv, sep="\t")

    kept_src = m[['kept_id', 'kept_source']].drop_duplicates()

    m_ctg = m[m['dup_source'] == 'contigs'].copy()
    agg_ctg = (
        m_ctg.groupby('kept_id')['dup_id']
            .agg(['count', lambda x: ",".join(x.astype(str).tolist()[:10])])
            .rename(columns={'count': 'n_trinity_dup', '<lambda_0>': 'trinity_list'})
            .reset_index()
    )
    agg_ctg['has_trinity_dup'] = agg_ctg['n_trinity_dup'] > 0

    m_ref = m[m['dup_source'] == 'ref'].copy()
    m_ref = m_ref[m_ref['dup_id'] != m_ref['kept_id']]
    agg_ref = (
        m_ref.groupby('kept_id')['dup_id']
            .agg(['count', lambda x: ",".join(x.astype(str).tolist()[:10])])
            .rename(columns={'count': 'n_ref_dup', '<lambda_0>': 'ref_dup_list'})
            .reset_index()
    )
    agg_ref['has_ref_dup'] = agg_ref['n_ref_dup'] > 0

    out = res_df.copy()
    out = out.merge(kept_src, left_on='Name', right_on='kept_id', how='left').drop(columns=['kept_id'])
    out = out.merge(agg_ctg,   left_on='Name', right_on='kept_id', how='left').drop(columns=['kept_id'])
    out = out.merge(agg_ref,   left_on='Name', right_on='kept_id', how='left').drop(columns=['kept_id'])

    # Normalize dtypes (avoid pandas FutureWarning)
    for col in ['kept_source', 'trinity_list', 'ref_dup_list']:
        if col in out.columns:
            out[col] = out[col].astype('string').fillna('')

    for col in ['n_trinity_dup', 'n_ref_dup']:
        if col in out.columns:
            out[col] = out[col].astype('Int64').fillna(0)

    for col in ['has_trinity_dup', 'has_ref_dup']:
        if col in out.columns:
            out[col] = out[col].astype('boolean').fillna(False)

    return out


def find_aberrant_noncoding(
    q_tumor: pd.DataFrame,
    q_ctrl: pd.DataFrame,
    known_noncoding_ids: Optional[Set[str]] = None,  # set of ENST IDs for known lncRNAs
    pseudo: float = 0.1,
    min_tpm_tumor: float = 1.0,
    max_tpm_ctrl: float = 0.5,
    min_log2fc: float = 2.0,
    top_n: Optional[int] = None,
    mark_trinity_novel_as_noncoding: bool = True
) -> pd.DataFrame:
    """
    Return columns:
      Name, TPM_tumor, TPM_ctrl, log2FC, is_aberrant, is_trinity, tri_ref_tx, enst_tag, nc_class
    Adds annotations only; does not filter types.
    """
    if known_noncoding_ids is None:
        known_noncoding_ids = set()

    t = q_tumor[['Name', 'TPM']].rename(columns={'TPM': 'TPM_tumor'})
    c = q_ctrl[['Name', 'TPM']].rename(columns={'TPM': 'TPM_ctrl'})
    df = t.merge(c, on='Name', how='inner')

    first_tok = df['Name'].str.split().str[0]
    is_trinity = first_tok.str.startswith('TRINITY')
    df['is_trinity'] = is_trinity

    tri_ref_tx = first_tok.where(is_trinity, None).apply(
        lambda x: (x.split('|')[1] if isinstance(x, str) and '|' in x else None)
    )
    df['tri_ref_tx'] = pd.Series(tri_ref_tx, index=df.index, dtype="object")

    def _pick_enst(tok, tri2):
        if isinstance(tok, str) and tok.startswith('TRINITY'):
            if isinstance(tri2, str) and tri2.startswith('ENST'):
                return tri2
            return None
        else:
            if isinstance(tok, str) and tok.startswith('ENST'):
                return tok
            return None

    df['enst_tag'] = [
        _pick_enst(ft, tr) for ft, tr in zip(first_tok.tolist(), df['tri_ref_tx'].tolist())
    ]

    enst_is_lnc = df['enst_tag'].isin(known_noncoding_ids)
    has_enst = df['enst_tag'].notna()

    nc_class = np.where(
        has_enst,
        np.where(enst_is_lnc, 'noncoding', 'coding_or_other'),
        np.where(
            df['tri_ref_tx'].eq('NOVEL') & mark_trinity_novel_as_noncoding,
            'novel',
            'unknown'
        )
    )
    df['nc_class'] = nc_class

    df['log2FC'] = np.log2((df['TPM_tumor'] + pseudo) / (df['TPM_ctrl'] + pseudo))
    df['is_aberrant'] = (
        (df['TPM_tumor'] >= min_tpm_tumor) &
        (df['TPM_ctrl'] <= max_tpm_ctrl) &
        (df['log2FC'] >= min_log2fc)
    )

    df = df.sort_values(['is_aberrant', 'log2FC', 'TPM_tumor'], ascending=[False, False, False])
    if top_n is not None:
        df = df.head(top_n)
    return df


def write_selected_pep(know_allow, novel_allow, pep_known, pep_novel, out_pep):
    """Write selected peptide FASTA entries (single-line sequences) to output."""
    out_pep = Path(out_pep)
    out_pep.parent.mkdir(parents=True, exist_ok=True)

    allow_known = set(know_allow)
    allow_novel = set(novel_allow)

    n_write = 0
    with open(out_pep, "w") as fo:
        for h, s in fasta_iter(pep_known):
            if _base_name_from_header(h) in allow_known:
                fo.write(">" + h + "\n")
                fo.write(s + "\n")  # write as a single line
                n_write += 1
        for h, s in fasta_iter(pep_novel):
            if _base_name_from_header(h) in allow_novel:
                fo.write(">" + h + "\n")
                fo.write(s + "\n")  # write as a single line
                n_write += 1
    print(f"[INFO] written {n_write} entries -> {out_pep}")


# -------------------------
# CLI
# -------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Select sORF-encoded peptides by aberrant noncoding expression (known + novel)."
    )
    # inputs
    ap.add_argument("--quant-tumor", required=True, help="salmon quant.sf for the tumor sample")
    ap.add_argument("--quant-ctrl",  required=True, help="salmon quant.sf for the control sample")
    ap.add_argument("--dedup-map",   required=True, help="merged_tx.dedup.map.tsv")
    ap.add_argument("--known-tx-fa", required=True, help="FASTA of known lncRNA transcripts (to load ENST set)")
    ap.add_argument("--pep-known",   required=True, help="FASTA of known sORF peptides")
    ap.add_argument("--pep-novel",   required=True, help="FASTA of novel sORF peptides")

    # thresholds
    ap.add_argument("--pseudo", type=float, default=0.1, help="pseudocount for log2FC calculation")
    ap.add_argument("--min-tpm-tumor", type=float, default=1.0, help="minimum tumor TPM")
    ap.add_argument("--max-tpm-ctrl",  type=float, default=0.5, help="maximum control TPM")
    ap.add_argument("--min-log2fc",    type=float, default=2.0, help="minimum log2 fold-change")
    ap.add_argument("--top-n", type=int, default=None, help="optional: keep top N rows after sorting")

    ap.add_argument(
        "--mark-trinity-novel-as-noncoding",
        action="store_true",
        default=True,
        help="treat TRINITY entries with NOVEL tag as noncoding (default: True)"
    )
    ap.add_argument(
        "--no-mark-trinity-novel-as-noncoding",
        dest="mark_trinity_novel_as_noncoding",
        action="store_false",
        help="disable treating TRINITY NOVEL as noncoding"
    )

    # outputs
    ap.add_argument("--out-dir",  required=True, help="output directory")
    ap.add_argument("--out-name", default="aeSEPs.pep", help="output peptide filename")
    ap.add_argument("--save-table", default=None, help="optional CSV path (relative to --out-dir) to save annotated table")

    return ap.parse_args()


def main():
    a = parse_args()

    # Input checks
    for p in [a.quant_tumor, a.quant_ctrl, a.dedup_map, a.known_tx_fa, a.pep_known, a.pep_novel]:
        check_file(p)

    # Load ENST set
    lnc_ids = load_fasta_ids(a.known_tx_fa)
    print(f"[INFO] loaded known lncRNA IDs: {len(lnc_ids)}")

    # Read quant tables
    q   = pd.read_csv(a.quant_tumor, sep="\t")
    q_c = pd.read_csv(a.quant_ctrl,  sep="\t")

    # Annotate aberrant/noncoding
    res = find_aberrant_noncoding(
        q, q_c,
        known_noncoding_ids=lnc_ids,
        pseudo=a.pseudo,
        min_tpm_tumor=a.min_tpm_tumor,
        max_tpm_ctrl=a.max_tpm_ctrl,
        min_log2fc=a.min_log2fc,
        top_n=a.top_n,
        mark_trinity_novel_as_noncoding=a.mark_trinity_novel_as_noncoding
    )
    res_anno = annotate_with_dedup_map(res, a.dedup_map)

    # Optional: save annotated table
    if a.save_table:
        Path(a.out_dir).mkdir(parents=True, exist_ok=True)
        out_csv = Path(a.out_dir) / a.save_table
        res_anno.to_csv(out_csv, index=False)
        print(f"[INFO] annotated table -> {out_csv}")

    # Filters for peptide selection
    know_res_anno  = res_anno[(res_anno['is_trinity'] == False) & (res_anno['is_aberrant'] == True) & (res_anno['nc_class'] == 'noncoding')]
    novel_res_anno = res_anno[(res_anno['is_trinity'] == True)  & (res_anno['is_aberrant'] == True)]

    # Select peptides
    out_pep = Path(a.out_dir) / a.out_name
    write_selected_pep(
        know_allow=know_res_anno["Name"].astype(str).tolist(),
        novel_allow=novel_res_anno["Name"].astype(str).tolist(),
        pep_known=a.pep_known,
        pep_novel=a.pep_novel,
        out_pep=out_pep
    )

    print("[DONE]")


if __name__ == "__main__":
    main()
