import os
import pandas as pd
from typing import Tuple, Iterable

# ---- Constants ----
# Default HLA class I gene set. You can override this by passing a custom set into `is_class_I`.
HLA_CLASS_I_GENES = {"A", "B", "C", "E", "F", "G"}


# ---- Utilities ----
def save_csv(MHCI: pd.DataFrame, MHCII: pd.DataFrame, save_file: str) -> Tuple[str, str]:
    """
    Save two DataFrames to CSV as <base>_MHCI.csv and <base>_MHCII.csv.

    Args:
        MHCI:  DataFrame for MHC class I results.
        MHCII: DataFrame for MHC class II results.
        save_file: Path ending with .csv; the suffix is used to derive the two output filenames.

    Returns:
        (mhci_csv_file, mhcii_csv_file): Paths to the saved CSV files.
    """
    if not save_file.lower().endswith(".csv"):
        raise ValueError("`save_file` must end with '.csv' to derive output names.")

    base_dir = os.path.dirname(os.path.abspath(save_file))
    if base_dir and not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    mhci_csv_file = save_file.replace(".csv", "_MHCI.csv")
    mhcii_csv_file = save_file.replace(".csv", "_MHCII.csv")

    MHCI.to_csv(mhci_csv_file, index=False)
    MHCII.to_csv(mhcii_csv_file, index=False)

    return mhci_csv_file, mhcii_csv_file


def is_class_I(allele: str, classI_genes: Iterable[str] = HLA_CLASS_I_GENES) -> bool:
    """
    Determine whether an allele belongs to MHC class I, based on its gene prefix.

    Rules:
      - Strip an optional "HLA-" prefix.
      - Split on '*' and take the gene part (e.g., "A*02:01" -> "A").
      - Check membership in the provided class I gene set.

    Args:
        allele: HLA allele string (e.g., "HLA-A*02:01", "A*02:01").
        classI_genes: Set-like of class I gene symbols.

    Returns:
        True if the allele maps to a class I gene; False otherwise.
    """
    if not isinstance(allele, str):
        return False
    a = allele.strip().replace("HLA-", "")
    gene = a.split("*", 1)[0]
    return gene in set(classI_genes)


def process_mhc_allele(allele_str: str, upper: bool = True) -> list[str]:
    """
    Normalize an allele string into a list of canonical "HLA-<gene>*<fields>" tokens.

    Example:
        "A*02:01/B*07:02" -> ["HLA-A*02:01", "HLA-B*07:02"]

    Args:
        allele_str: Raw allele string with '/' as delimiter.
        upper: If True, convert the entire token to uppercase.

    Returns:
        List of processed allele tokens.
    """
    if not allele_str:
        return []
    parts = [p.strip() for p in allele_str.split("/") if p.strip()]
    processed = []
    for p in parts:
        tok = p if p.startswith("HLA-") else f"HLA-{p}"
        processed.append(tok.upper() if upper else tok)
    return processed


def select_positive(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter rows whose 'Score' is not equal to the literal string 'Negative'.

    Note:
        - This is a simple string-based filter, matching your original logic.
        - If 'Score' can be numeric or categorical, adapt the predicate accordingly.

    Raises:
        KeyError if 'Score' column is missing.
    """
    if "Score" not in df.columns:
        raise KeyError("Column 'Score' is required for select_positive().")
    out = df.loc[df["Score"] != "Negative"].copy()
    out.reset_index(drop=True, inplace=True)
    return out


def sort_by_priority(df: pd.DataFrame, second_key: str) -> pd.DataFrame:
    """
    Sort by a predefined 'Score' priority (High > Moderate > Weak > Negative),
    then by a secondary key (ascending). Returns a copy without the helper column.

    Args:
        df: Input DataFrame; must contain 'Score' and `second_key`.
        second_key: Column name used as the tie-breaker (ascending).

    Returns:
        Sorted DataFrame (copy).

    Raises:
        KeyError if required columns are missing.
    """
    if "Score" not in df.columns:
        raise KeyError("Column 'Score' is required for sort_by_priority().")
    if second_key not in df.columns:
        raise KeyError(f"Column '{second_key}' is required for sort_by_priority().")

    priority_map = {"High": 3, "Moderate": 2, "Weak": 1, "Negative": 0}
    out = df.copy()
    out["_PriorityKey"] = out["Score"].map(priority_map).fillna(-1)

    out.sort_values(by=["_PriorityKey", second_key], ascending=[False, True], inplace=True)
    out.drop(columns=["_PriorityKey"], inplace=True)
    out.reset_index(drop=True, inplace=True)
    return out
