import os
import re
import pandas as pd


def filter_row(row: pd.Series, score_threshold: int) -> bool:
    """
    Filter alignment rows based on a simple priority rule:
      - Always keep rows already marked as best match.
      - Otherwise, parse the trailing integer score from `otherReferences` (comma-separated)
        and keep if score <= threshold.

    Args:
        row: Row from the alignment DataFrame. Must contain 'bestMatch' and 'otherReferences'.
        score_threshold: Minimum allowed alignment score threshold (lower is better).

    Returns:
        True if the row should be kept; False otherwise.
    """
    if bool(row.get("bestMatch", False)):
        return True  # keep rows explicitly marked as best match

    entry = str(row.get("otherReferences", ""))
    parts = entry.split(",")
    if len(parts) >= 4:
        score_str = parts[-1].strip()
        try:
            score = int(score_str)
            return score <= score_threshold
        except ValueError:
            # Malformed score field; drop conservatively
            print("[filter_row] warning: invalid score field ->", score_str)
            return False
    else:
        # Malformed 'otherReferences' structure; drop conservatively
        print("[filter_row] warning: malformed otherReferences ->", entry)
        return False


def sum_numbers_before_M(s: str) -> int:
    """
    Sum all numbers that precede 'M' in a CIGAR-like string.

    Example:
        "10M2I5M" -> 10 + 5 = 15

    Args:
        s: CIGAR segment string.

    Returns:
        Sum of all integer counts directly preceding 'M'.
    """
    numbers = re.findall(r'(\d+)M', str(s))
    return sum(int(num) for num in numbers)


def contains_match_length(row: pd.Series,
                          best_match_length: int,
                          match_length_threshold: float) -> bool:
    """
    Check if the last-but-one comma field of 'otherReferences' contains a CIGAR-like
    segment whose total 'M' length is >= best_match_length * match_length_threshold.

    Args:
        row: Row containing 'otherReferences'.
        best_match_length: Length of the best (primary) match (sum of 'M').
        match_length_threshold: Fractional threshold to accept alternative matches.

    Returns:
        True if the alternative match satisfies the length requirement; otherwise False.
    """
    other = str(row.get('otherReferences', ''))
    if 'M' not in other:
        return False
    # Expecting "...,<CIGAR>,<score>"; we grab the CIGAR part with [-2]
    parts = other.split(',')
    if len(parts) < 2:
        return False
    cigar = parts[-2]
    m_len = sum_numbers_before_M(cigar)
    return m_len >= best_match_length * match_length_threshold


def exploded_other_references(df: pd.DataFrame,
                              best_match_length: int,
                              match_length_threshold: float,
                              pair: bool) -> pd.DataFrame:
    """
    Expand semicolon-separated 'otherReferences' into individual rows,
    retain those with sufficient match length, and append back to the original
    best-match rows.

    Args:
        df: Input alignment DataFrame; must contain 'otherReferences' and core columns.
        best_match_length: Length of best (primary) match.
        match_length_threshold: Fractional threshold to accept secondary matches.
        pair: Whether data are paired-end (kept for interface parity; not used here).

    Returns:
        Combined DataFrame with original best matches and accepted alternative matches.
    """
    base = df.copy()
    base['bestMatch'] = True

    exploded = base.copy()
    exploded['otherReferences'] = exploded['otherReferences'].astype(str).str.split(';')
    exploded = exploded.explode('otherReferences')
    # Keep only the portion after the last ':' to unify format
    exploded['otherReferences'] = exploded['otherReferences'].astype(str).str.split(':').str[-1]

    mask = exploded.apply(
        contains_match_length,
        axis=1,
        args=(best_match_length, match_length_threshold)
    )
    matched_df = exploded.loc[mask].copy()

    # Parse fields from "nucleotideId,matchPos,matchLength,score" format
    parts = matched_df['otherReferences'].str.split(',')
    matched_df['nucleotideId'] = parts.str[0]
    matched_df['matchPos'] = parts.str[1]
    matched_df['matchLength'] = parts.str[2]
    matched_df['numbers'] = matched_df['matchLength'].apply(sum_numbers_before_M)
    matched_df['bestMatch'] = False

    combined = pd.concat([base, matched_df], ignore_index=True)
    return combined


def get_match_df(sample: str,
                 output_nucleic: str,
                 match_length_threshold: float,
                 pair: bool) -> tuple[pd.DataFrame, int]:
    """
    Parse a per-sample alignment summary (<sample>.txt), compute the best-match length,
    and keep rows whose 'M' sum is within the threshold of the best.

    Expected input columns (tab-delimited):
        [0] read_id
        [2] reference
        [3] position
        [5] CIGAR-like length string
        [6] next_reference
        [9] read_sequence
        [11] other_reference blob

    Args:
        sample: Sample identifier.
        output_nucleic: Directory containing "<sample>.txt".
        match_length_threshold: Fractional threshold.
        pair: Paired-end indicator (not used here; kept for interface parity).

    Returns:
        (final_match_df, best_match_length)
    """
    reads_ids, reads = [], []
    references, next_references = [], []
    match_length, match_pos, other_references = [], [], []

    txt_file = os.path.join(output_nucleic, f"{sample}.txt")
    with open(txt_file, 'r') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if not line:
                continue
            cols = line.split("\t")
            # guard against short lines
            if len(cols) < 12:
                continue

            read_id = cols[0]
            reference = cols[2]
            pos = cols[3]
            length = cols[5]
            next_reference = cols[6]
            read_seq = cols[9]
            other_reference = cols[11]

            if reference == '*':
                continue

            reads_ids.append(read_id)
            reads.append(read_seq)
            references.append(reference)
            next_references.append(next_reference)
            other_references.append(other_reference)
            match_length.append(length)
            match_pos.append(pos)

    match_df = pd.DataFrame(
        {
            'readsId': reads_ids,
            'reads': reads,
            'nucleotideId': references,
            'nextReferences': next_references,
            'matchLength': match_length,
            'matchPos': match_pos,
            'otherReferences': other_references,
        }
    )
    match_df['numbers'] = match_df['matchLength'].apply(sum_numbers_before_M)

    best_match_length = int(match_df['numbers'].max()) if not match_df.empty else 0
    final_match_df = match_df.loc[
        match_df['numbers'] >= best_match_length * match_length_threshold
    ].copy()
    final_match_df['bestMatch'] = True
    return final_match_df, best_match_length


def get_microbe_abundance(muilt_match_df_merge: pd.DataFrame,
                          tax_id_hierarchy_df: pd.DataFrame,
                          total_library_size: int) -> pd.DataFrame:
    """
    Compute per-tax_id abundance metrics:
      - Reads: number of unique (readsId, reads) pairs
      - Best Match Reads: count of rows flagged as bestMatch
      - Score: sum over 1/matchTimes
      - eRPKM: Score normalized by genome length (kb) and library size (M)

    Args:
        muilt_match_df_merge: Merged alignment + taxonomy + genome length table.
        tax_id_hierarchy_df: Taxonomy hierarchy table.
        total_library_size: Total library size (host + non-host reads).

    Returns:
        DataFrame with abundance and taxonomy columns, plus eRPKM Normalized per superkingdom.
    """
    micro_abundance_list = []
    for tax_id, dfg in muilt_match_df_merge.groupby('tax_id'):
        # Unique reads per tax_id
        reads_cnt = dfg[['readsId', 'reads']].drop_duplicates().shape[0]
        best_match_reads = dfg.loc[dfg['bestMatch'] == True].shape[0]
        score_sum = (1 / dfg['matchTimes']).sum()
        genome_length = dfg['genomeLength'].iloc[0]
        erpkm = score_sum / (genome_length / 1e3) / (total_library_size / 1e6)
        micro_abundance_list.append([
            tax_id, genome_length, total_library_size, reads_cnt,
            best_match_reads, score_sum, erpkm
        ])

    micro_abundance = pd.DataFrame(
        micro_abundance_list,
        columns=[
            'Tax ID', 'Genome Length', 'Total Library Size',
            'Reads', 'Best Match Reads', 'Score', 'eRPKM'
        ]
    )

    # Compute unique/ambiguous read counts per tax_id
    rmdup = muilt_match_df_merge.drop_duplicates(
        subset=['readsId', 'reads', 'tax_id'], keep='first'
    ).copy()
    match_times = rmdup.groupby(['readsId', 'reads']).size().reset_index(name='matchTimes')
    rmdup = rmdup.merge(match_times, on=['readsId', 'reads'], how='left')

    unambiguous = (
        rmdup.assign(_is_unique=lambda x: (x['matchTimes'] == 1).astype(int))
             .groupby('tax_id')['_is_unique'].sum()
             .reset_index()
             .rename(columns={'_is_unique': 'Unique Match Reads', 'tax_id': 'Tax ID'})
    )

    # Merge taxonomy and unique-read stats
    micro_abundance = micro_abundance.merge(
        tax_id_hierarchy_df, left_on='Tax ID', right_on='tax_id', how='left'
    ).merge(unambiguous, on='Tax ID', how='left')

    final_header = [
        'Tax ID', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom',
        'tax_id_scientific_name', 'species_scientific_name', 'genus_scientific_name',
        'family_scientific_name', 'order_scientific_name', 'class_scientific_name',
        'phylum_scientific_name', 'superkingdom_scientific_name',
        'Genome Length', 'Total Library Size', 'Reads', 'Unique Match Reads',
        'Best Match Reads', 'Score', 'eRPKM'
    ]
    micro_abundance = micro_abundance.loc[:, final_header]

    # Normalize eRPKM within each superkingdom (percentage)
    micro_abundance['eRPKM Normalized'] = (
        micro_abundance.groupby("superkingdom_scientific_name", group_keys=False)["eRPKM"]
        .transform(lambda x: (x / x.sum()) * 100 if x.sum() else 0)
    )
    return micro_abundance


def get_microbe_rank_abundance(muilt_match_df_merge: pd.DataFrame,
                               microbe_abundance_df: pd.DataFrame,
                               rank_name: str) -> pd.DataFrame:
    """
    Aggregate abundance metrics at a specified taxonomy rank (e.g., 'species', 'genus').

    Args:
        muilt_match_df_merge: Merged alignment table.
        microbe_abundance_df: Base abundance table returned by `get_microbe_abundance`.
        rank_name: Taxonomic rank name to aggregate by (column present in taxonomy table).

    Returns:
        Rank-level abundance DataFrame with normalized eRPKM.
    """
    rank_cols = [f'{rank_name}', f"{rank_name}_scientific_name", "superkingdom_scientific_name"]
    sum_cols = ['Best Match Reads', 'Score', 'eRPKM']
    non_sum_cols = ['Total Library Size']

    rank_df = microbe_abundance_df[rank_cols + sum_cols + non_sum_cols].copy()
    rank_df = rank_df.groupby(rank_cols + non_sum_cols, as_index=False).sum()

    # Unique/ambiguous counts per read at the given rank
    rmdup = muilt_match_df_merge.drop_duplicates(
        subset=['readsId', 'reads', f'{rank_name}'], keep='first'
    ).copy()
    match_times = rmdup.groupby(['readsId', 'reads']).size().reset_index(name='matchTimes')
    rmdup = rmdup.merge(match_times, on=['readsId', 'reads'], how='left')

    reads_list = []
    for tax_id, dfg in rmdup.groupby(f'{rank_name}'):
        reads_cnt = dfg.shape[0]
        unique_cnt = (dfg['matchTimes'] == 1).sum()
        reads_list.append([tax_id, reads_cnt, unique_cnt])

    reads_df = pd.DataFrame(reads_list, columns=[f'{rank_name}', 'Reads', 'Unique Match Reads'])
    rank_df = rank_df.merge(reads_df, how='left', on=f'{rank_name}')

    # Normalize eRPKM within superkingdom at this rank
    rank_df['eRPKM Normalized'] = (
        rank_df.groupby("superkingdom_scientific_name", group_keys=False)["eRPKM"]
        .transform(lambda x: (x / x.sum()) * 100 if x.sum() else 0)
    )

    # Cast rank ID to int if possible (guarded)
    try:
        rank_df[rank_name] = rank_df[rank_name].astype(int)
    except Exception:
        pass

    final_header = [
        f"{rank_name}", f"{rank_name}_scientific_name", "superkingdom_scientific_name",
        "Total Library Size", "Reads", "Unique Match Reads",
        "Best Match Reads", "Score", "eRPKM", "eRPKM Normalized"
    ]
    rank_df = rank_df.loc[:, final_header]
    rank_df = rank_df.rename(columns={
        f"{rank_name}": f"{rank_name.capitalize()} Tax Id",
        f"{rank_name}_scientific_name": f"{rank_name.capitalize()} Scientific Name",
        "superkingdom_scientific_name": "Superkingdom Scientific Name"
    })
    return rank_df


def get_tax_id(sample: str,
               match_df: pd.DataFrame,
               catalog_genome_rank_df: pd.DataFrame,
               tax_id_hierarchy_df: pd.DataFrame,
               output_nucleic: str) -> pd.DataFrame:
    """
    Merge matches to genome rank catalog and taxonomy, extract unique tax_ids,
    and write `<sample>.txids`.

    Args:
        sample: Sample ID.
        match_df: Match table from `get_match_df`.
        catalog_genome_rank_df: Catalog with 'genome_accession' and 'tax_id'.
        tax_id_hierarchy_df: Hierarchy table with 'tax_id' and 'species_scientific_name'.
        output_nucleic: Output directory.

    Returns:
        The merged DataFrame with taxonomy columns added.
    """
    merged = match_df.merge(
        catalog_genome_rank_df, left_on='nucleotideId', right_on='genome_accession', how='left'
    ).merge(
        tax_id_hierarchy_df[['tax_id', 'species_scientific_name']],
        on='tax_id', how='left'
    )

    unique_ids = merged['tax_id'].dropna().drop_duplicates()
    txids_path = os.path.join(output_nucleic, f"{sample}.txids")
    unique_ids.to_csv(txids_path, header=None, index=None)
    return merged


def get_fa(sample: str,
           final_match_df_merge: pd.DataFrame,
           output_nucleic: str) -> None:
    """
    Write a FASTA with headers carrying mapping metadata.

    Header format:
        >ntIndex:{i}|<readsId>|<matchLength>|<matchPos>|<ref>|<taxId>|<speciesName_underscored>

    Args:
        sample: Sample ID.
        final_match_df_merge: Merged table containing required columns.
        output_nucleic: Output directory.
    """
    lines = []
    n = final_match_df_merge.shape[0]
    for i in range(n):
        reads = final_match_df_merge['reads'].iloc[i]
        readsId = final_match_df_merge['readsId'].iloc[i]
        matchLength = final_match_df_merge['matchLength'].iloc[i]
        matchPos = str(final_match_df_merge['matchPos'].iloc[i])
        ref = final_match_df_merge['nucleotideId'].iloc[i]
        taxId = str(final_match_df_merge['tax_id'].iloc[i])
        speciesName = str(final_match_df_merge['species_scientific_name'].iloc[i]).replace(" ", "_")
        header = f">ntIndex:{i}|{readsId}|{matchLength}|{matchPos}|{ref}|{taxId}|{speciesName}"
        lines.append(header)
        lines.append(reads)

    sseq_file = os.path.join(output_nucleic, f"{sample}.fasta")
    with open(sseq_file, "w") as f:
        f.write("\n".join(lines))


def microbe_abundancing(sample: str,
                        match_df: pd.DataFrame,
                        output_nucleic: str,
                        best_match_length: int,
                        match_length_threshold: float,
                        score_threshold: int,
                        catalog_genome_rank_df: pd.DataFrame,
                        tax_id_hierarchy_df: pd.DataFrame,
                        microbial_genome_length_with_taxid_sum_df: pd.DataFrame,
                        total_library_size: int,
                        pair: bool) -> None:
    """
    Full abundance workflow:
      1) Expand 'otherReferences' into alternative matches and filter by length + score.
      2) Deduplicate read-reference pairs and compute per-read matchTimes.
      3) Merge with genome/taxonomy/length tables.
      4) Compute abundance at tax_id and rank levels.
      5) Persist intermediate and final CSVs.

    Args: see parameter names; all required.
    """
    # 1) explode and filter by length/score
    muilt_match_df = exploded_other_references(
        match_df, best_match_length, match_length_threshold, pair
    )
    muilt_match_df_filtered = muilt_match_df[
        muilt_match_df.apply(filter_row, axis=1, args=(score_threshold,))
    ].copy()

    # 2) deduplicate and compute matchTimes
    muilt_match_df_filtered.drop_duplicates(
        subset=['readsId', 'reads', 'nucleotideId'], keep='first', inplace=True
    )
    match_times = muilt_match_df_filtered.groupby(['readsId', 'reads']).size().reset_index(name='matchTimes')
    muilt_match_df_filtered = muilt_match_df_filtered.merge(match_times, on=['readsId', 'reads'])

    # 3) taxonomy / genome merges
    muilt_match_df_merge = muilt_match_df_filtered.merge(
        catalog_genome_rank_df, left_on='nucleotideId', right_on='genome_accession', how='left'
    )
    muilt_match_df_merge.pop('genome_accession')
    muilt_match_df_merge = muilt_match_df_merge.merge(
        microbial_genome_length_with_taxid_sum_df, how='left'
    )
    muilt_match_df_merge = muilt_match_df_merge.merge(
        tax_id_hierarchy_df, left_on='tax_id', right_on='tax_id', how='left'
    )

    # 4) abundance summaries
    microbe_abundance_df = get_microbe_abundance(
        muilt_match_df_merge, tax_id_hierarchy_df, total_library_size
    )
    microbe_abundance_df_species = get_microbe_rank_abundance(
        muilt_match_df_merge, microbe_abundance_df, 'species'
    )
    microbe_abundance_df_genus = get_microbe_rank_abundance(
        muilt_match_df_merge, microbe_abundance_df, 'genus'
    )

    # 5) write outputs
    muilt_match_df_merge.to_csv(os.path.join(output_nucleic, f"{sample}.muilt_match.csv"), index=False)
    microbe_abundance_df.to_csv(os.path.join(output_nucleic, f"{sample}.microbe_abundance.csv"), index=False)
    microbe_abundance_df_species.to_csv(os.path.join(output_nucleic, f"{sample}.microbe_abundance_species.csv"), index=False)
    microbe_abundance_df_genus.to_csv(os.path.join(output_nucleic, f"{sample}.microbe_abundance_genus.csv"), index=False)


def get_data(sample: str,
             catalog_genome_rank_file: str,
             tax_id_hierarchy_file: str,
             output_hg38: str,
             output_nucleic: str,
             microbial_genome_length: str,
             match_length_threshold: float,
             score_threshold: int,
             pair: bool,
             seq_type: str) -> None:
    """
    Main I/O wrapper that drives the per-sample analysis:
      - Read taxonomy and genome length catalogs
      - Compute total library size from flagstat
      - Parse matches, attach taxonomy, export FASTA
      - Run abundance summaries and write CSVs

    Args:
        sample: Sample identifier.
        catalog_genome_rank_file: TSV with columns including 'genome_accession' and 'tax_id'.
        tax_id_hierarchy_file: TSV taxonomy hierarchy.
        output_hg38: Directory containing `<sample>_<seq_type>_hg38_align_sort.bam.flagstat.txt`.
        output_nucleic: Directory for microbial outputs.
        microbial_genome_length: CSV with columns: nucleotideId,genomeLength (no header variants).
        match_length_threshold: Ratio threshold for alternative match acceptance.
        score_threshold: Score cutoff for `filter_row`.
        pair: Paired-end indicator.
        seq_type: Sequencing type label used in file naming.
    """
    catalog_genome_rank_df = pd.read_csv(catalog_genome_rank_file, sep="\t")
    tax_id_hierarchy_df = pd.read_csv(tax_id_hierarchy_file, sep="\t")

    # Load genome lengths and sum by tax_id
    gl_df = pd.read_csv(microbial_genome_length, sep=',')
    gl_df.columns = ['nucleotideId', 'genomeLength']
    gl_tax = gl_df.merge(
        catalog_genome_rank_df, left_on='nucleotideId', right_on='genome_accession', how='left'
    ).drop(columns=['nucleotideId', 'genome_accession'])
    microbial_genome_length_with_taxid_sum_df = (
        gl_tax.groupby("tax_id", as_index=False)['genomeLength'].sum()
              .rename(columns={'genomeLength': 'genomeLength'})
    )

    # Total library size from samtools flagstat file
    flagstat_path = os.path.join(output_hg38, f"{sample}_{seq_type}_hg38_align_sort.bam.flagstat.txt")
    with open(flagstat_path, 'r') as f:
        first_line = f.readline().strip()
        total_library_size = int(first_line.split(' ')[0])

    # Parse matches and attach taxonomy
    match_df, best_match_length = get_match_df(sample, output_nucleic, match_length_threshold, pair)
    match_df_merge = get_tax_id(sample, match_df, catalog_genome_rank_df, tax_id_hierarchy_df, output_nucleic)

    # Export FASTA with metadata headers
    get_fa(sample, match_df_merge, output_nucleic)

    # Abundance summaries
    microbe_abundancing(
        sample=sample,
        match_df=match_df,
        output_nucleic=output_nucleic,
        best_match_length=best_match_length,
        match_length_threshold=match_length_threshold,
        score_threshold=score_threshold,
        catalog_genome_rank_df=catalog_genome_rank_df,
        tax_id_hierarchy_df=tax_id_hierarchy_df,
        microbial_genome_length_with_taxid_sum_df=microbial_genome_length_with_taxid_sum_df,
        total_library_size=total_library_size,
        pair=pair
    )
