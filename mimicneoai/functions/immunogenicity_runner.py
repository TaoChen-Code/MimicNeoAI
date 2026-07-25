from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from mimicneoai.immunogenicity_prediction.core import InferenceConfig, run_inference
from mimicneoai.immunogenicity_prediction.runtime.defaults import run_default_inference



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


def predict_default_immunogenicity_df(
    df: pd.DataFrame,
    *,
    antigen_class: str,
    hla_source: str = "netmhc-pseudoseq",
    hla_fasta_path: str = "",
    hla_pseudoseq_csv: Sequence[str] = (),
    peptide_col: str = "peptide",
    hla_col: str = "hla",
    output_score_col: str = "immunogenicity_score",
    batch_size: int = 512,
    device: str = "auto",
    num_processes: int = 1,
    include_input_qc: bool = True,
    model_root: str | None = None,
) -> pd.DataFrame:
    """Run a recovered default model selected by antigen class.

    The cryptic class is evaluated as the packaged ensemble mean. The default
    HLA source matches the recovered runtime checkpoint metadata.
    """
    cfg = InferenceConfig(
        model_path="",
        hla_fasta_path=hla_fasta_path,
        hla_source=hla_source,
        hla_pseudoseq_csv=tuple(hla_pseudoseq_csv),
        peptide_col=peptide_col,
        hla_col=hla_col,
        output_score_col=output_score_col,
        batch_size=batch_size,
        device=device,
        num_processes=num_processes,
        include_input_qc=include_input_qc,
    )
    return run_default_inference(df, antigen_class, cfg, model_root=model_root)


def predict_default_immunogenicity_csv(
    input_csv: str,
    output_csv: str,
    *,
    antigen_class: str,
    hla_source: str = "netmhc-pseudoseq",
    hla_fasta_path: str = "",
    hla_pseudoseq_csv: Sequence[str] = (),
    peptide_col: str = "peptide",
    hla_col: str = "hla",
    output_score_col: str = "immunogenicity_score",
    batch_size: int = 512,
    device: str = "auto",
    num_processes: int = 1,
    include_input_qc: bool = True,
    model_root: str | None = None,
) -> pd.DataFrame:
    """CSV counterpart of :func:`predict_default_immunogenicity_df`."""
    output = predict_default_immunogenicity_df(
        pd.read_csv(input_csv),
        antigen_class=antigen_class,
        hla_source=hla_source,
        hla_fasta_path=hla_fasta_path,
        hla_pseudoseq_csv=hla_pseudoseq_csv,
        peptide_col=peptide_col,
        hla_col=hla_col,
        output_score_col=output_score_col,
        batch_size=batch_size,
        device=device,
        num_processes=num_processes,
        include_input_qc=include_input_qc,
        model_root=model_root,
    )
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_csv, index=False)
    return output
