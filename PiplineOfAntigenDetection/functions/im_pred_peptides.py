# import
import os
import pandas as pd
import warnings
from multiprocessing import Pool, cpu_count
warnings.filterwarnings("ignore", category=FutureWarning)  # 屏蔽所有 FutureWarning

# calculate_peptide_features

def calculate_peptide_features(df: pd.DataFrame, peptide_col: str, R_HOME, R_LIBRARY) -> pd.DataFrame:
    """
    在 Python 中调用 R 的 Peptides 包计算肽段理化特征
    
    参数：
        df : 包含肽段序列的 Pandas 数据框
        peptide_col : 数据框中用于存储肽段序列的列名
        
    返回：
        添加了 8 个新特征列的数据框：
        aa_composition, polarity, volume, net_charge, 
        hydrophobicity, boman_index, aliphatic_index, isoelectric_point
    """    

    # 设置 R 的环境变量
    os.environ['R_HOME'] = R_HOME

    # 导入 R 相关库
    import anndata2ri
    import rpy2.robjects
    import rpy2.robjects as ro
    from rpy2.robjects import r, pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
    
    # 激活 pandas 到 R 数据框的转换
    pandas2ri.activate()
    anndata2ri.activate()

    # 引入 R 的 Peptides 包
    # 设置自定义 R 库路径
    r(f'.libPaths(c("{R_LIBRARY}"))')
    
    # 获取 R 环境
    base = importr('base')
    utils = importr('utils')
    Peptides = importr('Peptides')
    
    # 在 R 中定义函数
    r_code = f'''
    calculate_features <- function(data, peptide_col) {{
        suppressPackageStartupMessages(library(Peptides))
        
        calculate_polarity <- function(seq) {{
            polarities <- list(
                hydrophilic = c('A','D','E','N','Q','R','S','T','Y'),
                hydrophobic = c('I','L','M','F','P','V','W')
            )
            sum(sapply(strsplit(seq, '')[[1]], function(aa) {{
                if (aa %in% polarities$hydrophilic) 1
                else if (aa %in% polarities$hydrophobic) -1
                else 0
            }}))
        }}
        
        aa_volumes <- c(
            A=88.6, C=118.8, D=111.1, E=138.3, F=189.9, G=60.1,
            H=153.2, I=166.6, K=168.6, L=166.6, M=162.9, N=114.1,
            P=115.0, Q=146.2, R=174.0, S=105.9, T=119.0, V=140.0,
            W=227.8, Y=193.6
        )
        
        calculate_volume <- function(seq) {{
            sum(sapply(strsplit(seq, '')[[1]], function(aa) {{
                if (aa %in% names(aa_volumes)) aa_volumes[aa] else 0
            }}))
        }}
        
        data$aa_composition <- sapply(data[[peptide_col]], function(x) 
            paste(unlist(aaComp(x)), collapse=", "))
        data$polarity <- sapply(data[[peptide_col]], calculate_polarity)
        data$volume <- sapply(data[[peptide_col]], calculate_volume)
        data$net_charge <- sapply(data[[peptide_col]], charge)
        data$hydrophobicity <- sapply(data[[peptide_col]], hydrophobicity)
        data$boman_index <- sapply(data[[peptide_col]], boman)
        data$aliphatic_index <- sapply(data[[peptide_col]], aIndex)
        data$isoelectric_point <- sapply(data[[peptide_col]], pI)
        
        return(data)
    }}
    '''
    
    # 执行 R 代码，定义函数
    ro.r(r_code)
    
    # 获取 R 中的 calculate_features 函数
    calculate_features = ro.globalenv['calculate_features']
    
    # 转换并执行计算
    with (ro.default_converter + pandas2ri.converter).context():
        r_df = ro.conversion.get_conversion().py2rpy(df)
        result_r = calculate_features(r_df, peptide_col)
        result_df = ro.conversion.get_conversion().rpy2py(result_r)
    
    return result_df

def parallel_calculate_features(R_HOME, R_LIBRARY, peptides, seq_col, num_chunks=1):
    """
    确保所有行被分配，最后一个 chunk 包含剩余行
    
    Args:
        peptides (pd.DataFrame): 输入的肽段数据
        seq_col (str): 肽序列列名
        num_chunks (int): 拆分的 chunk 数量（也是进程数）
        
    Returns:
        pd.DataFrame: 合并后的特征计算结果
    """
    # 计算每个 chunk 的基础大小和余数
    total_rows = len(peptides)
    chunk_size = total_rows // num_chunks
    remainder = total_rows % num_chunks
    
    # 拆分 chunks
    chunks = []
    start = 0
    for i in range(num_chunks):
        # 前 'remainder' 个 chunk 多分配 1 行（均匀分配余数）
        end = start + chunk_size + (1 if i < remainder else 0)
        chunks.append(peptides.iloc[start:end])
        start = end
    
    # 进程数 = chunk 数
    with Pool(processes=len(chunks)) as pool:
        results = pool.starmap(
            calculate_peptide_features,
            [(chunk, seq_col, R_HOME, R_LIBRARY) for chunk in chunks]
        )
    ###
    return pd.concat(results, ignore_index=True)