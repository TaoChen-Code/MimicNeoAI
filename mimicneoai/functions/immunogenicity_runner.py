from __future__ import annotations

from typing import Optional

import pandas as pd

from mimicneoai.immunogenicity_prediction.core import InferenceConfig, run_inference



def predict_immunogenicity_df(
    df: pd.DataFrame,
    *,
    model_path: str,
    hla_fasta_path: str,
    peptide_col: str = "peptide",
    hla_col: str = "hla",
    output_score_col: str = "immunogenicity_score",
    batch_size: int = 512,
    device: str = "auto",
    num_processes: int = 1,
) -> pd.DataFrame:
    cfg = InferenceConfig(
        model_path=model_path,
        hla_fasta_path=hla_fasta_path,
        peptide_col=peptide_col,
        hla_col=hla_col,
        output_score_col=output_score_col,
        batch_size=batch_size,
        device=device,
        num_processes=num_processes,
    )
    return run_inference(df, cfg)



def predict_immunogenicity_csv(
    input_csv: str,
    output_csv: str,
    *,
    model_path: str,
    hla_fasta_path: str,
    peptide_col: str = "peptide",
    hla_col: str = "hla",
    output_score_col: str = "immunogenicity_score",
    batch_size: int = 512,
    device: str = "auto",
    num_processes: int = 1,
) -> pd.DataFrame:
    df = pd.read_csv(input_csv)
    out = predict_immunogenicity_df(
        df,
        model_path=model_path,
        hla_fasta_path=hla_fasta_path,
        peptide_col=peptide_col,
        hla_col=hla_col,
        output_score_col=output_score_col,
        batch_size=batch_size,
        device=device,
        num_processes=num_processes,
    )
    out.to_csv(output_csv, index=False)
    return out
