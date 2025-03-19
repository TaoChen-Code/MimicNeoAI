import os
import pandas as pd

def erpkm(df,sample,type):
    # Standardize eRPKM values and rename columns
    df = df.rename(columns = {"eRPKM Normalized":sample,f'{type.capitalize()} Scientific Name':'name'})
    return df[['name',sample]].reset_index(drop=True)

def get_data_for_stack_plot(df,sampleName,topN):
    # Convert wide-format dataframe to long-format for visualization
    df_filtered = df.head(topN)
    category = [sampleName for _ in range(df_filtered.shape[0]+1)]
    subcategory = df_filtered['name'].tolist() + ['others']
    value = df_filtered[sampleName].tolist()
    value.append(100-sum(value))
    return category, subcategory, value

def merge_and_mean(samples,seq_type,type,group):
    # Aggregate and average data across multiple samples
    df = pd.DataFrame()
    for sample in samples:
        temp_df = pd.read_csv(f"{ipt_dir}/{sample}/04.nucleic_v12/{sample}.microbe_abundance_{type}.csv")
        temp_df = erpkm(temp_df,sample,type)
        df = pd.merge(df, temp_df, how='outer') if not df.empty else temp_df
    
    df = df.fillna(0)
    df[group] = df.iloc[:,1:].mean(axis=1)
    df[group] = df[group] / df[group].sum() * 100
    return df[['name',group]].sort_values(group, ascending=False)

def transfer(group1_samples,group2_samples,group1_name,group2_name,seq_type,type,topN):
    # Process data for comparative analysis between two groups
    group1 = merge_and_mean(group1_samples,seq_type,type,group1_name)
    group2 = merge_and_mean(group2_samples,seq_type,type,group2_name)
    
    merged = pd.merge(group1, group2, on='name', how='outer').fillna(0)
    merged = merged.sort_values(group1_name, ascending=False)
    
    # Prepare visualization data
    g1_cat, g1_subcat, g1_val = get_data_for_stack_plot(merged[[group1_name]], group1_name, topN)
    g2_cat, g2_subcat, g2_val = get_data_for_stack_plot(merged[[group2_name]], group2_name, topN)
    
    return pd.DataFrame({
        'category': g1_cat + g2_cat,
        'subcategory': g1_subcat + g2_subcat,
        'value': g1_val + g2_val
    })

# Response group samples (renamed)
group1_samples = ["P0001", "P0004", "P0006", "P0010"]
# Non-response group samples (renamed)
group2_samples = ["P0002", "P0003", "P0005", "P0009"]
ipt_dir = "./data/"
seq_type = 'rna'

# Generate genus-level comparison
genus_comparison = transfer(group1_samples, group2_samples, 
                           "Response", "Non_Response", seq_type, "genus", 10)
# Generate species-level comparison
species_comparison = transfer(group1_samples, group2_samples, 
                             "Response", "Non_Response", seq_type, "species", 10)

# Export results
genus_comparison.to_csv("./results/Response_Non-Response_genus.txt", sep="\t", index=False)
species_comparison.to_csv("./results/Response_Non-Response_species.txt", sep="\t", index=False)