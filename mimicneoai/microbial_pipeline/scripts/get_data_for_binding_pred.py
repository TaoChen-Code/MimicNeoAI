import os
import pandas as pd

def parse_yp_taxids(qseqid: str):
    """
    Parse tax_ids from YP:Z: field in qseqid.
    Example qseqid:
      'B1|E2500...|99|150M|YP:Z:1912894,28116,371601,1912896'
    Returns:
      ['1912894', '28116', '371601', '1912896']
    """
    if not isinstance(qseqid, str):
        return []
    if "YP:Z:" not in qseqid:
        return []
    tail = qseqid.split("YP:Z:", 1)[1]
    tail = tail.split("|", 1)[0]
    return [x.strip() for x in tail.split(",") if x.strip()]


def filter_by_catalog_blastx(df_blast: pd.DataFrame, catalog_df: pd.DataFrame) -> pd.DataFrame:
    """
    Catalog filter for BLASTX results, analogous to DIAMOND filtering:
      - Normalize BLASTX sseqid (e.g. 'ref|WP_231809274.1|') to 'WP_231809274.1'
      - Parse YP tax_ids from qseqid
      - Expand to (qseqid, sseqid_clean, yp_taxid)
      - Inner-join with catalog (prot_id, tax_id)
      - Collapse back to rows that have at least one matched tax_id
      - Append matched tax_ids to qseqid as '|CTAXA:tax1,tax2,...'

    Returns a filtered DataFrame.
    """
    df = df_blast.copy()

    # catalog: prot_id, tax_id (tax_id as str)
    catalog_sub = catalog_df[["prot_id", "tax_id"]].copy()
    catalog_sub["tax_id"] = catalog_sub["tax_id"].astype(str)

    # normalize BLASTX sseqid, e.g. 'ref|WP_231809274.1|' -> 'WP_231809274.1'
    df["sseqid_clean"] = (
        df["sseqid"]
        .astype(str)
        .str.replace(r"^ref\|", "", regex=True)  # drop leading 'ref|'
        .str.rstrip("|")                         # drop trailing '|'
    )

    # parse YP taxids from qseqid
    df["yp_taxid_list"] = df["qseqid"].apply(parse_yp_taxids)

    # expand to long format: (qseqid, sseqid_clean, yp_taxid)
    tmp = df[["qseqid", "sseqid_clean", "yp_taxid_list"]].explode("yp_taxid_list")
    tmp = tmp.rename(columns={"yp_taxid_list": "yp_taxid"})
    tmp = tmp.dropna(subset=["yp_taxid"])

    if tmp.empty:
        return df.iloc[0:0].copy()

    # inner-join with catalog on (sseqid_clean, yp_taxid) ↔ (prot_id, tax_id)
    merged = pd.merge(
        tmp,
        catalog_sub,
        left_on=["sseqid_clean", "yp_taxid"],
        right_on=["prot_id", "tax_id"],
        how="inner",
    )
    if merged.empty:
        return df.iloc[0:0].copy()

    # collect matched tax_ids per (qseqid, sseqid_clean)
    matched = (
        merged.groupby(["qseqid", "sseqid_clean"])["yp_taxid"]
        .agg(lambda xs: sorted(set(xs)))
        .reset_index()
        .rename(columns={"yp_taxid": "matched_taxids"})
    )

    # merge back to original df (still keeps original sseqid for context)
    df_merged = df.merge(matched, on=["qseqid", "sseqid_clean"], how="left")
    df_filtered = df_merged[df_merged["matched_taxids"].notna()].copy()
    if df_filtered.empty:
        return df_filtered

    # build CTAXA suffix and append to qseqid
    df_filtered["matched_taxids_str"] = df_filtered["matched_taxids"].apply(
        lambda lst: ",".join(lst) if isinstance(lst, list) else str(lst)
    )

    # build CTAXA suffix and append to qseqid (add protein id too)
    df_filtered["qseqid"] = (
            df_filtered["qseqid"]
            + "|CTAXA:" + df_filtered["matched_taxids_str"]
            + "|PROT:" + df_filtered["sseqid_clean"].astype(str)
    )

    # drop intermediate columns
    df_filtered = df_filtered.drop(
        columns=["yp_taxid_list", "matched_taxids", "matched_taxids_str", "sseqid_clean"],
        errors="ignore",
    )

    return df_filtered



def get_data_for_binding_pred(
    blast_file: str,
    colnames: list,
    pvacbind_file: str,
    output_blastx: str,
    min_pident: float,
    catalog_df: pd.DataFrame,
    max_evalue: float = 1e-5,
    min_qcovs: float = 90.0,   # <-- renamed
    sample=None,
    tool=None,
):
    """
    Process BLASTX results and generate peptide FASTA for binding prediction (pVACbind).

    Steps:
      1) Read BLASTX TSV file.
      2) Apply catalog filter (PathSeq YP tax_ids must be supported by V7_catalog_prot4).
      3) Filter by pident >= min_pident.
      4) Filter by evalue <= max_evalue.
      5) Filter by BLASTX qcovs >= min_qcovs (official subject coverage metric).
      6) Extract final peptides (sseq) into FASTA.
      7) Write summary table.

    Args:
        blast_file (str): BLASTX result filename.
        colnames (list): Column names from outfmt 6.
        pvacbind_file (str): Output FASTA file.
        output_blastx (str): Directory of BLASTX outputs.
        min_pident (float): Minimum identity.
        catalog_df (DataFrame): V7_catalog_prot4.
        max_evalue (float): E-value threshold.
        min_qcovs (float): Minimum qcovs threshold.
        sample (str|None): Sample ID.
        tool: Logger with write_log.

    Returns:
        dict: Summary statistics.
    """
    if tool is None:
        class DummyTool:
            def write_log(self, msg, level="info"):
                print(f"[{level.upper()}] {msg}")
        tool = DummyTool()

    blast_path = os.path.join(output_blastx, blast_file)
    if not os.path.exists(blast_path):
        tool.write_log(f"[{sample}] BLAST file not found: {blast_path}", "warning")
        return {}

    # Read BLASTX table
    df = pd.read_csv(blast_path, sep="\t", header=None, names=colnames)

    summary = {}

    def _log_step(name, d):
        summary[name] = len(d)
        tool.write_log(f"[{sample}] {name}: {len(d)} hits", "info")

    # raw
    _log_step("raw_blastx", df)

    if df.empty:
        open(os.path.join(output_blastx, pvacbind_file), "w").close()
        return summary

    # convert numeric columns
    for col in ("pident", "evalue", "qcovs"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # step: after_taxonomic_consistency_check
    # BLASTX protein hits are retained only if their microbial taxon IDs
    # match PathSeq-derived YP taxon IDs at the species level, ensuring
    # species-level taxonomic consistency between nucleic and protein evidence.
    df_cat = filter_by_catalog_blastx(df, catalog_df)
    _log_step("after_taxonomic_consistency_check", df_cat)
    if df_cat.empty:
        open(os.path.join(output_blastx, pvacbind_file), "w").close()
        return summary

    # pident filtering
    df_pid = df_cat[df_cat["pident"] >= float(min_pident)].copy()
    _log_step(f"after_pident_{int(min_pident)}", df_pid)
    if df_pid.empty:
        open(os.path.join(output_blastx, pvacbind_file), "w").close()
        return summary

    # evalue filtering
    df_eval = df_pid[df_pid["evalue"] <= float(max_evalue)].copy()
    _log_step(f"after_evalue_{max_evalue}", df_eval)
    if df_eval.empty:
        open(os.path.join(output_blastx, pvacbind_file), "w").close()
        return summary

    # qcovs filtering
    if "qcovs" in df_eval.columns:
        df_qcov = df_eval[df_eval["qcovs"] >= float(min_qcovs)].copy()
    else:
        tool.write_log(
            f"[{sample}] WARNING: qcovs column missing; skipping coverage filter",
            "warning",
        )
        df_qcov = df_eval.copy()

    _log_step(f"after_qcovs_{int(min_qcovs)}", df_qcov)

    # write summary
    summary_path = os.path.join(
        output_blastx,
        f"{sample if sample else 'blastx'}.binding_pred_filter_stats.tsv"
    )
    pd.DataFrame(
        [{"step": k, "n_hits": v} for k, v in summary.items()]
    ).to_csv(summary_path, sep="\t", index=False)

    if df_qcov.empty:
        open(os.path.join(output_blastx, pvacbind_file), "w").close()
        return summary

    # ------------------------------------------------------------
    # Assign protIndex and write filtered BLASTX table + FASTA
    # ------------------------------------------------------------
    # 1) Create a copy to avoid side-effects
    df_qcov_out = df_qcov.copy().reset_index(drop=True)

    # 2) Add protIndex column (1-based) and rewrite qseqid as:
    #    protIndex:{i}|{original_qseqid_with_Ctaxa_and_Prot}
    df_qcov_out["protIndex"] = range(1, len(df_qcov_out) + 1)
    df_qcov_out["qseqid"] = df_qcov_out.apply(
        lambda r: f"protIndex:{int(r['protIndex'])}|{r['qseqid']}",
        axis=1
    )

    # 3) Save filtered BLASTX table for downstream merging
    filtered_path = os.path.join(output_blastx, f"{sample}.blastx.filtered")
    df_qcov_out.to_csv(filtered_path, sep="\t", index=False)
    tool.write_log(f"[{sample}] Wrote filtered BLASTX table: {filtered_path}", "info")

    # 4) Write FASTA using the SAME qseqid (already prefixed with protIndex)
    fasta_lines = []
    for row in df_qcov_out.itertuples(index=False):
        header = f">{row.qseqid}"
        peptide = str(row.sseq).replace("*", "")
        if not peptide:
            continue
        fasta_lines.append(header)
        fasta_lines.append(peptide)

    fasta_path = os.path.join(output_blastx, pvacbind_file)
    with open(fasta_path, "w") as f_out:
        f_out.write("\n".join(fasta_lines))

    tool.write_log(
        f"[{sample}] Wrote {len(fasta_lines) // 2} peptides to {fasta_path}",
        "info"
    )

    return summary

