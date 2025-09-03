import os
import pandas as pd

# tools
def save_csv(MHCI, MHCII, save_file):
    mhci_csv_file = save_file.replace('.csv', '_MHCI.csv')  # 生成 MHCI CSV 文件名
    mhcii_csv_file = save_file.replace('.csv', '_MHCII.csv')  # 生成 MHCII CSV 文件名
    
    MHCI.to_csv(mhci_csv_file, index=False)
    MHCII.to_csv(mhcii_csv_file, index=False)

    return mhci_csv_file,mhcii_csv_file

def is_class_I(allele):
    # 去掉前缀 "HLA-"
    allele = allele.replace("HLA-", "")
    # 以星号分割，获得主要基因部分
    gene = allele.split('*')[0]

    # 判断是否属于一类
    return gene in hlaI_set

def process_mhc_allele(allele_str):
    alleles = allele_str.split("/")
    processed = []
    for allele in alleles:
        if not allele.startswith("HLA-"):
            allele = f"HLA-{allele}"
        processed.append(allele)
    return processed

def select_positive(df):
    return df[df['Score'] != 'Negative'].reset_index(drop=True)

def sort_by_priority(df,second_key):
    # 创建优先级映射字典
    priority = {'High': 3, 'Moderate': 2, 'Weak': 1,'Negative':0}
    df['Priority'] = df['Score'].map(priority)
    
    # 按Epitope Seq分组，保留每组中优先级最高的记录
    sorted_df = df.sort_values(by=['Priority',second_key], ascending=[False,True])
    sorted_df.pop('Priority')
    return sorted_df