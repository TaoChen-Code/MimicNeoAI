import pandas as pd
import numpy as np
import os
import re

def filter_row(row, score_threshold):
    """Filter alignment rows based on match quality
    
    Args:
        row (pd.Series): Row from alignment DataFrame
        score_threshold (int): Minimum allowed alignment score threshold
        
    Returns:
        bool: True if row meets filtering criteria, False otherwise"""
    if row["bestMatch"]:
        return True  # Keep rows marked as best match
    
    entry = str(row["otherReferences"])
    parts = entry.split(",")
    if len(parts) >= 4: 
        score_str = parts[-1].strip()  
        try:
            score = int(score_str)
            if score <= score_threshold:
                return True
        except ValueError:
            print("filter_row error!")
            return False
    else:
        print("filter_row error!")
    return False

def sum_numbers_before_M(s):
    """Calculate sum of numbers preceding 'M' in CIGAR string
    
    Args:
        s (str): CIGAR format string containing matches
        
    Returns:
        int: Sum of all numbers preceding 'M' characters"""
    numbers = re.findall(r'(\d*)M',s)
    return sum(int(num) for num in numbers)

def contains_match_length(row, best_match_length, match_length_threshold):
    """Check if alternative references contain sufficient match length
    
    Args:
        row (pd.Series): Row from alignment DataFrame
        best_match_length (int): Length of best match
        match_length_threshold (float): Ratio threshold for acceptable match length
        
    Returns:
        bool: True if secondary matches meet length requirements"""
    if 'M' not in str(row['otherReferences']):
        return False
    M_match_length = sum_numbers_before_M(row['otherReferences'].split(',')[-2])
    return M_match_length >= best_match_length * match_length_threshold

def exploded_other_references(df, best_match_length, match_length_threshold, pair):
    """Expand and filter alternative alignment references
    
    Args:
        df (pd.DataFrame): Input alignment data
        best_match_length (int): Length of best match sequence
        match_length_threshold (float): Minimum match length ratio threshold
        pair (bool): Indicates paired-end sequencing data
        
    Returns:
        pd.DataFrame: Expanded DataFrame with valid alternative references"""
    df['bestMatch'] = True
    df['otherReferences'] = df['otherReferences'].str.split(';')
    df_exploded = df.explode('otherReferences')
    df_exploded['otherReferences'] = df_exploded['otherReferences'].astype(str)
    df_exploded['otherReferences'] = df_exploded['otherReferences'].str.split(':').str[-1]
    
    mask = df_exploded.apply(contains_match_length,axis=1, args=(best_match_length,match_length_threshold,))
    matched_df = df_exploded[mask].copy()
    matched_df['nucleotideId'] = matched_df['otherReferences'].str.split(',').str[0]
    matched_df['matchPos'] = matched_df['otherReferences'].str.split(',').str[1]
    matched_df['matchLength'] = matched_df['otherReferences'].str.split(',').str[2]
    matched_df['numbers'] = matched_df['matchLength'].apply(sum_numbers_before_M)
    matched_df['bestMatch'] = False
    
    df_combined = pd.concat([df, matched_df], ignore_index=True)
    
    return df_combined

def get_match_df(sample, output_nucleic, match_length_threshold, pair):
    """Process alignment file and extract match information
    
    Args:
        sample (str): Sample identifier
        output_nucleic (str): Output directory path
        match_length_threshold (float): Minimum match length ratio threshold
        pair (bool): Indicates paired-end sequencing data
        
    Returns:
        tuple: (Filtered matches DataFrame, best match length)"""
    references,reads_ids,reads,match_length,next_references,other_references,match_pos = [],[],[],[],[],[],[]
    with open(f"{output_nucleic}/{sample}.txt", 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line_split = line.split("\t")
        id = line_split[0]
        reference = line_split[2]
        pos = line_split[3]
        length = line_split[5]
        next_reference = line_split[6]
        read = line_split[9]
        other_reference = line_split[11]
        if reference == '*':
            continue
        reads_ids.append(id)
        reads.append(read)
        references.append(reference)
        next_references.append(next_reference)
        other_references.append(other_reference)
        match_length.append(length)
        match_pos.append(pos)
    
    match_df = pd.DataFrame(zip(reads_ids,reads,references,next_references,match_length,match_pos,other_references),columns = ['readsId','reads','nucleotideId','nextReferences','matchLength','matchPos','otherReferences'])
    match_df['numbers'] = match_df['matchLength'].apply(sum_numbers_before_M)
    best_match_length = max(match_df['numbers'])
    final_match_df = match_df[match_df['numbers'] >= best_match_length * match_length_threshold].copy()
    final_match_df['bestMatch'] = True
    return final_match_df, best_match_length

def get_microbe_abundance(muilt_match_df_merge, tax_id_hierarchy_df, total_library_size):
    """Calculate microbial abundance metrics
    
    Args:
        muilt_match_df_merge (pd.DataFrame): Merged alignment data
        tax_id_hierarchy_df (pd.DataFrame): Taxonomic hierarchy data
        total_library_size (int): Total sequencing library size
        
    Returns:
        pd.DataFrame: Abundance metrics with taxonomic information"""
    micro_abundance_list = []
    for (tax_id,df) in muilt_match_df_merge.groupby('tax_id'):        
        reads = df[['readsId', 'reads']].drop_duplicates().shape[0]
        best_match_reads = df[df['bestMatch'] == True].shape[0]
        muilt_match_reads_sum = (1 / df['matchTimes']).sum()
        genome_length = df['genomeLength'].iloc[0]
        erpkm = muilt_match_reads_sum / (genome_length/1e3) / (total_library_size/1e6)
        micro_abundance_list.append([tax_id,genome_length,total_library_size,reads,best_match_reads,muilt_match_reads_sum,erpkm])
    micro_abundance = pd.DataFrame(micro_abundance_list,columns=['Tax ID','Genome Length','Total Library Size','Reads','Best Match Reads','Score','eRPKM'])
    
    muilt_match_df_merge_rmdup = muilt_match_df_merge.drop_duplicates(subset=['readsId','reads','tax_id'], keep='first').copy()
    match_times = muilt_match_df_merge_rmdup.groupby(['readsId','reads']).size().reset_index(name='matchTimes')
    muilt_match_df_merge_rmdup.pop('matchTimes')
    muilt_match_df_merge_rmdup_merged = muilt_match_df_merge_rmdup.merge(match_times, on=['readsId','reads'])
    
    unambiguous_list = []
    for (tax_id,df) in muilt_match_df_merge_rmdup_merged.groupby('tax_id'):
        unambiguous_reads = df[df['matchTimes'] == 1].shape[0]
        unambiguous_list.append([tax_id,unambiguous_reads]) 
    unambiguous_df = pd.DataFrame(unambiguous_list,columns=['Tax ID','Unique Match Reads'])
    
    micro_abundance = pd.merge(micro_abundance,tax_id_hierarchy_df,left_on='Tax ID',right_on='tax_id',how='left')
    micro_abundance = pd.merge(micro_abundance,unambiguous_df,on='Tax ID',how='left')
    final_header = ['Tax ID','species','genus','family','order','class','phylum','superkingdom','tax_id_scientific_name',"species_scientific_name",'genus_scientific_name',
                    'family_scientific_name','order_scientific_name','class_scientific_name','phylum_scientific_name','superkingdom_scientific_name',
                    'Genome Length','Total Library Size','Reads','Unique Match Reads','Best Match Reads','Score','eRPKM']
    micro_abundance = micro_abundance.loc[:,final_header]
    micro_abundance['eRPKM Normalized'] = (
        micro_abundance.groupby("superkingdom_scientific_name", group_keys=False)["eRPKM"]
        .transform(lambda x: x / x.sum() * 100)
    )
    return micro_abundance

def get_microbe_rank_abundance(muilt_match_df_merge, microbe_abundance_df, rank_name):
    """Calculate abundance metrics for specified taxonomic rank
    
    Args:
        muilt_match_df_merge (pd.DataFrame): Merged alignment data
        microbe_abundance_df (pd.DataFrame): Base abundance metrics
        rank_name (str): Taxonomic rank to aggregate (e.g., 'species')
        
    Returns:
        pd.DataFrame: Aggregated abundance metrics for specified rank"""
    rank_id_and_name = [f'{rank_name}',f"{rank_name}_scientific_name","superkingdom_scientific_name"]
    sum_info_columns = ['Best Match Reads','Score','eRPKM']
    non_sum_info_columns = ['Total Library Size']
    microbe_abundance_df_rank = microbe_abundance_df[rank_id_and_name + sum_info_columns + non_sum_info_columns].copy()

    microbe_abundance_df_rank = microbe_abundance_df_rank.groupby(rank_id_and_name + non_sum_info_columns).sum()
    microbe_abundance_df_rank = microbe_abundance_df_rank.reset_index()
    
    muilt_match_df_merge_rmdup = muilt_match_df_merge.drop_duplicates(subset=['readsId','reads',f'{rank_name}'], keep='first').copy()
    match_times = muilt_match_df_merge_rmdup.groupby(['readsId','reads']).size().reset_index(name='matchTimes')
    muilt_match_df_merge_rmdup.pop('matchTimes')
    muilt_match_df_merge_rmdup_merged = muilt_match_df_merge_rmdup.merge(match_times, on=['readsId','reads'])
    
    reads_list = []
    for (tax_id,df) in muilt_match_df_merge_rmdup_merged.groupby(f'{rank_name}'):
        reads = df.shape[0]
        unambiguous_reads = df[df['matchTimes'] == 1].shape[0]
        reads_list.append([tax_id,reads,unambiguous_reads])
    reads_df = pd.DataFrame(reads_list,columns=[f'{rank_name}','Reads','Unique Match Reads'])
    microbe_abundance_df_rank = pd.merge(microbe_abundance_df_rank,reads_df,how='left')
    
    microbe_abundance_df_rank['eRPKM Normalized'] = (
        microbe_abundance_df_rank.groupby("superkingdom_scientific_name", group_keys=False)["eRPKM"]
        .transform(lambda x: x / x.sum() * 100)
    )
    microbe_abundance_df_rank[rank_name] = microbe_abundance_df_rank[rank_name].astype(int)

    final_header = [f"{rank_name}",f"{rank_name}_scientific_name","superkingdom_scientific_name","Total Library Size","Reads","Unique Match Reads",
                    "Best Match Reads","Score","eRPKM","eRPKM Normalized"]					
    microbe_abundance_df_rank = microbe_abundance_df_rank.loc[:,final_header]
    microbe_abundance_df_rank = microbe_abundance_df_rank.rename(columns={
                                                                 f"{rank_name}":f"{rank_name.capitalize()} Tax Id",
                                                                 f"{rank_name}_scientific_name":f"{rank_name.capitalize()} Scientific Name",
                                                                 "superkingdom_scientific_name":"Superkingdom Scientific Name"})
    return microbe_abundance_df_rank

def get_tax_id(sample, match_df, catalog_genome_rank_df, tax_id_hierarchy_df, output_nucleic):
    """Extract taxonomic IDs and write to file
    
    Args:
        sample (str): Sample identifier
        match_df (pd.DataFrame): Alignment matches data
        catalog_genome_rank_df (pd.DataFrame): Genome taxonomy catalog
        tax_id_hierarchy_df (pd.DataFrame): Taxonomic hierarchy data
        output_nucleic (str): Output directory path"""
    match_df_merge = pd.merge(match_df,catalog_genome_rank_df,left_on='nucleotideId',right_on='genome_accession',how='left')
    match_df_merge = pd.merge(match_df_merge,tax_id_hierarchy_df[['tax_id','species_scientific_name']],on='tax_id',how='left')
    unique_match_df = match_df_merge.drop_duplicates(subset='tax_id')
    unique_match_df['tax_id'].to_csv(f"{output_nucleic}/{sample}.txids",header=None,index=None)
    return match_df_merge

def get_fa(sample, final_match_df_merge, output_nucleic):
    """Generate FASTA file from alignment results
    
    Args:
        sample (str): Sample identifier
        final_match_df_merge (pd.DataFrame): Processed alignment data
        output_nucleic (str): Output directory path"""
    results = ''
    for i in range(final_match_df_merge.shape[0]):
        reads = final_match_df_merge['reads'].iloc[i]
        readsId = final_match_df_merge['readsId'].iloc[i]
        matchLength = final_match_df_merge['matchLength'].iloc[i]
        matchPos = str(final_match_df_merge['matchPos'].iloc[i])
        ref = final_match_df_merge['nucleotideId'].iloc[i]
        taxId = str(final_match_df_merge['tax_id'].iloc[i])
        speciesName = final_match_df_merge['species_scientific_name'].iloc[i].replace(" ","_")
        results += f">ntIndex:{i}|" + readsId +"|" + matchLength +"|" + matchPos + "|" + ref + "|" + taxId + "|" + speciesName + '\n'
        results += reads + '\n'
    results = results[:-1]
    sseq_file = f"{output_nucleic}/{sample}.fasta"
    with open(sseq_file,"w")as f:
        f.write(results)

def microbe_abundancing(sample, match_df, output_nucleic, best_match_length, match_length_threshold, 
                       score_threshold, catalog_genome_rank_df, tax_id_hierarchy_df, 
                       microbial_genome_length_with_taxid_sum_df, total_library_size, pair):
    """Main microbial abundance analysis workflow
    
    Args:
        sample (str): Sample identifier
        match_df (pd.DataFrame): Initial alignment matches
        output_nucleic (str): Output directory path
        best_match_length (int): Length of best alignment match
        match_length_threshold (float): Match length ratio threshold
        score_threshold (int): Alignment score cutoff
        catalog_genome_rank_df (pd.DataFrame): Genome taxonomy catalog
        tax_id_hierarchy_df (pd.DataFrame): Taxonomic hierarchy data
        microbial_genome_length_with_taxid_sum_df (pd.DataFrame): Genome length data
        total_library_size (int): Total sequencing depth (host + non-host reads).
        pair (bool): Paired-end data indicator"""
    
    muilt_match_df = exploded_other_references(match_df,best_match_length,match_length_threshold,pair)
    muilt_match_df_filtered = muilt_match_df[muilt_match_df.apply(filter_row, axis=1, args=(score_threshold,))]
    
    muilt_match_df_filtered = muilt_match_df_filtered.drop_duplicates(subset=['readsId','reads', 'nucleotideId'], keep='first')
    match_times = muilt_match_df_filtered.groupby(['readsId','reads']).size().reset_index(name='matchTimes')
    muilt_match_df_filtered = muilt_match_df_filtered.merge(match_times, on=['readsId','reads'])
    
    muilt_match_df_merge = pd.merge(muilt_match_df_filtered,catalog_genome_rank_df,left_on='nucleotideId',right_on='genome_accession',how='left')
    muilt_match_df_merge.pop('genome_accession')
    muilt_match_df_merge = pd.merge(muilt_match_df_merge,microbial_genome_length_with_taxid_sum_df,how='left')
    muilt_match_df_merge = pd.merge(muilt_match_df_merge,tax_id_hierarchy_df,left_on='tax_id',right_on='tax_id',how='left')
    
    microbe_abundance_df = get_microbe_abundance(muilt_match_df_merge,tax_id_hierarchy_df,total_library_size)
    microbe_abundance_df_species = get_microbe_rank_abundance(muilt_match_df_merge,microbe_abundance_df,'species')
    microbe_abundance_df_genus = get_microbe_rank_abundance(muilt_match_df_merge,microbe_abundance_df,'genus')
    
    muilt_match_df_merge.to_csv(f"{output_nucleic}/{sample}.muilt_match.csv",index=None)
    microbe_abundance_df.to_csv(f"{output_nucleic}/{sample}.microbe_abundance.csv",index=None)
    microbe_abundance_df_species.to_csv(f"{output_nucleic}/{sample}.microbe_abundance_species.csv",index=None)
    microbe_abundance_df_genus.to_csv(f"{output_nucleic}/{sample}.microbe_abundance_genus.csv",index=None)

def get_data(sample, catalog_genome_rank_file, tax_id_hierarchy_file, output_hg38, output_nucleic,
            microbial_genome_length, match_length_threshold, score_threshold, pair, seq_type):
    """Main data processing pipeline
    
    Args:
        sample (str): Sample identifier
        catalog_genome_rank_file (str): Path to genome taxonomy catalog
        tax_id_hierarchy_file (str): Path to taxonomic hierarchy file
        output_hg38 (str): Human genome alignment output directory
        output_nucleic (str): Microbial analysis output directory
        microbial_genome_length (str): Path to microbial genome lengths file
        match_length_threshold (float): Match length ratio threshold
        score_threshold (int): Alignment score cutoff
        pair (bool): Paired-end data indicator
        seq_type (str): Sequencing type identifier"""
    catalog_genome_rank_df = pd.read_csv(catalog_genome_rank_file,sep="\t")
    tax_id_hierarchy_df = pd.read_csv(tax_id_hierarchy_file,sep="\t")
    microbial_genome_length_df = pd.read_csv(microbial_genome_length,sep=',')
    microbial_genome_length_df.columns=['nucleotideId','genomeLength']
    microbial_genome_length_with_taxid_df = pd.merge(microbial_genome_length_df,catalog_genome_rank_df,left_on='nucleotideId',right_on='genome_accession')
    microbial_genome_length_with_taxid_df.pop('nucleotideId')
    microbial_genome_length_with_taxid_df.pop('genome_accession')
    microbial_genome_length_with_taxid_sum_df = microbial_genome_length_with_taxid_df.groupby("tax_id").sum()
    microbial_genome_length_with_taxid_sum_df = microbial_genome_length_with_taxid_sum_df.reset_index()
    microbial_genome_length_with_taxid_sum_df = microbial_genome_length_with_taxid_sum_df.rename({'genomeLength':'Genome Length'})
    total_library = f"{output_hg38}/{sample}_{seq_type}_hg38_align_sort.bam.flagstat.txt"
    with open(total_library,'r') as f:
        contexts = f.readlines()
        total_library_size = int(contexts[0].split(' ')[0])
    
    match_df,best_match_length = get_match_df(sample,output_nucleic,match_length_threshold,pair)
    match_df_merge = get_tax_id(sample,match_df,catalog_genome_rank_df,tax_id_hierarchy_df,output_nucleic)
    get_fa(sample,match_df_merge,output_nucleic)
    microbe_abundancing(sample,match_df,output_nucleic,best_match_length,match_length_threshold,score_threshold,catalog_genome_rank_df,\
                        tax_id_hierarchy_df,microbial_genome_length_with_taxid_sum_df,total_library_size,pair)
