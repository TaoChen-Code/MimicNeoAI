# coding=utf-8
# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import os
import pandas as pd
import warnings
from multiprocessing import Pool, cpu_count

# Silence noisy FutureWarnings that are not actionable for callers
warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------
# calculate_peptide_features
# ---------------------------------------------------------------------
def calculate_peptide_features(
    df: pd.DataFrame,
    peptide_col: str,
    R_HOME: str,
    R_LIBRARY: str,
) -> pd.DataFrame:
    """
    Compute peptide physicochemical features via R's Peptides package (through rpy2).

    Args:
        df (pd.DataFrame): Input DataFrame containing peptide sequences.
        peptide_col (str): Column name in `df` that stores peptide sequences (AAs).
        R_HOME (str): Path to the R installation (used to set env var R_HOME).
        R_LIBRARY (str): Path to custom R library directory containing 'Peptides', etc.

    Returns:
        pd.DataFrame: Same rows as input with 8 additional feature columns:
            - aa_composition
            - polarity
            - volume
            - net_charge
            - hydrophobicity
            - boman_index
            - aliphatic_index
            - isoelectric_point
    """
    if df is None or len(df) == 0:
        return df

    if peptide_col not in df.columns:
        raise KeyError(f"Column `{peptide_col}` not found in DataFrame.")

    # Configure R environment
    if R_HOME:
        os.environ["R_HOME"] = R_HOME

    # Lazy import rpy2 stack (faster module import when not needed)
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects import r
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    # Enable pandas <-> R DataFrame automatic conversion
    pandas2ri.activate()
    anndata2ri.activate()

    # Set custom R library path (prepend to existing)
    r(f'.libPaths(c("{R_LIBRARY}", .libPaths()))')

    # Load required R packages
    base = importr("base")
    utils = importr("utils")
    Peptides = importr("Peptides")

    # R helper: feature calculation on a data.frame
    r_code = r"""
    calculate_features <- function(data, peptide_col) {
      suppressPackageStartupMessages(library(Peptides))

      # Simple polarity score (+1 hydrophilic, -1 hydrophobic, 0 otherwise)
      calculate_polarity <- function(seq) {
        polarities <- list(
          hydrophilic = c('A','D','E','N','Q','R','S','T','Y'),
          hydrophobic = c('I','L','M','F','P','V','W')
        )
        sum(sapply(strsplit(seq, '')[[1]], function(aa) {
          if (aa %in% polarities$hydrophilic) 1
          else if (aa %in% polarities$hydrophobic) -1
          else 0
        }))
      }

      # Per-residue volumes (A^3), coarse scale
      aa_volumes <- c(
        A=88.6, C=118.8, D=111.1, E=138.3, F=189.9, G=60.1,
        H=153.2, I=166.6, K=168.6, L=166.6, M=162.9, N=114.1,
        P=115.0, Q=146.2, R=174.0, S=105.9, T=119.0, V=140.0,
        W=227.8, Y=193.6
      )

      calculate_volume <- function(seq) {
        sum(sapply(strsplit(seq, '')[[1]], function(aa) {
          if (aa %in% names(aa_volumes)) aa_volumes[aa] else 0
        }))
      }

      # Feature columns (vectorized via sapply over the peptide column)
      data$aa_composition    <- sapply(data[[peptide_col]], function(x) paste(unlist(aaComp(x)), collapse=", "))
      data$polarity          <- sapply(data[[peptide_col]], calculate_polarity)
      data$volume            <- sapply(data[[peptide_col]], calculate_volume)
      data$net_charge        <- sapply(data[[peptide_col]], charge)
      data$hydrophobicity    <- sapply(data[[peptide_col]], hydrophobicity)
      data$boman_index       <- sapply(data[[peptide_col]], boman)
      data$aliphatic_index   <- sapply(data[[peptide_col]], aIndex)
      data$isoelectric_point <- sapply(data[[peptide_col]], pI)

      return(data)
    }
    """
    ro.r(r_code)
    calculate_features = ro.globalenv["calculate_features"]

    # Convert to R, compute, convert back
    with (ro.default_converter + pandas2ri.converter).context():
        r_df = ro.conversion.get_conversion().py2rpy(df)
        result_r = calculate_features(r_df, peptide_col)
        result_df = ro.conversion.get_conversion().rpy2py(result_r)

    return result_df


# ---------------------------------------------------------------------
# parallel_calculate_features
# ---------------------------------------------------------------------
def parallel_calculate_features(
    R_HOME: str,
    R_LIBRARY: str,
    peptides: pd.DataFrame,
    seq_col: str,
    num_chunks: int = 1,
) -> pd.DataFrame:
    """
    Split the input DataFrame into `num_chunks` (nearly equal-sized) parts,
    compute features in parallel, and concatenate the results.

    - All rows are assigned; the first `remainder` chunks receive one extra row.
    - Each worker calls `calculate_peptide_features` independently.

    Args:
        R_HOME (str): Path to the R installation (env var R_HOME).
        R_LIBRARY (str): Path to custom R library directory.
        peptides (pd.DataFrame): Input peptide table.
        seq_col (str): Column name containing peptide sequences.
        num_chunks (int): Number of chunks / processes to use.

    Returns:
        pd.DataFrame: Concatenated feature table across all chunks.
    """
    if peptides is None or len(peptides) == 0:
        return peptides

    if seq_col not in peptides.columns:
        raise KeyError(f"Column `{seq_col}` not found in DataFrame.")

    # Normalize chunk count
    if not isinstance(num_chunks, int) or num_chunks < 1:
        num_chunks = 1

    # Avoid oversubscription on tiny inputs
    num_chunks = min(num_chunks, len(peptides), cpu_count() or 1)

    total_rows = len(peptides)
    base = total_rows // num_chunks
    remainder = total_rows % num_chunks

    # Build chunks with even remainder distribution
    chunks = []
    start = 0
    for i in range(num_chunks):
        end = start + base + (1 if i < remainder else 0)
        chunks.append(peptides.iloc[start:end])
        start = end

    # Parallel execution (process-per-chunk)
    with Pool(processes=len(chunks)) as pool:
        results = pool.starmap(
            calculate_peptide_features,
            [(chunk, seq_col, R_HOME, R_LIBRARY) for chunk in chunks],
        )

    return pd.concat(results, ignore_index=True)
