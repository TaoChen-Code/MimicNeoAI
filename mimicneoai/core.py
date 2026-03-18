import re
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from tqdm import tqdm

from mimicneoai.functions.nodemon_pool import NoDaemonPool
from mimicneoai.immunogenicity_prediction.model import Vocab, build_model, load_model_weights


AMINO = "ARNDCQEGHILKMFPSTWYV-"
PADDING_LENGTH = 419
X2_DIM = 25

R_FEATURE_CODE = """
calculate_features <- function(data, peptide_col) {
    suppressPackageStartupMessages(library(Peptides))

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
}
"""


@dataclass
class InferenceConfig:
    model_path: str
    hla_fasta_path: str
    peptide_col: str = "peptide"
    hla_col: str = "hla"
    output_score_col: str = "immunogenicity_score"
    batch_size: int = 512
    device: str = "auto"
    num_processes: int = 1
    verbose: bool = True


def pick_device(device_name: str) -> torch.device:
    name = str(device_name).strip().lower()
    if name == "cpu":
        return torch.device("cpu")
    if name == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but not available.")
        return torch.device("cuda:0")
    if name.startswith("cuda:"):
        if not torch.cuda.is_available():
            raise RuntimeError(f"{device_name} requested but CUDA is not available.")
        idx = int(name.split(":", 1)[1])
        if idx < 0 or idx >= torch.cuda.device_count():
            raise RuntimeError(
                f"CUDA device index out of range: {idx}. "
                f"Available count: {torch.cuda.device_count()}."
            )
        return torch.device(f"cuda:{idx}")
    return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


def normalize_hla(hla: str) -> str:
    val = str(hla).strip()
    val = re.sub(r"^HLA-", "", val, flags=re.IGNORECASE)
    if "-" in val:
        val = re.sub(r"(?<=[0-9])-(?=[A-Z])", "/", val)
    if "/" in val:
        val = val.split("/")[-1]
    parts = val.split(":")
    if len(parts) >= 2:
        val = f"{parts[0]}:{parts[1]}"
    return val


def parse_hla_fasta(hla_fasta_path: str) -> Dict[str, str]:
    hla_names: List[str] = []
    prots: List[str] = []

    with open(hla_fasta_path, "r", encoding="utf-8") as f:
        prot = ""
        started = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                tokens = line.split()
                chosen = None
                for tok in tokens[1:]:
                    clean = tok.replace("HLA-", "")
                    if "*" in clean:
                        chosen = clean
                        break
                if chosen is None:
                    chosen = tokens[0].replace(">", "")
                hla_names.append(chosen)
                if started:
                    prots.append(prot)
                    prot = ""
                started = True
            else:
                prot += line

    if started:
        prots.append(prot)

    hla_df = pd.DataFrame({"HLA": hla_names, "protin": prots})
    hla_df["HLA"] = hla_df["HLA"].astype(str).str.split(":").str[:2].str.join(":")
    hla_df = hla_df.drop_duplicates(subset="HLA", keep="first").reset_index(drop=True)
    return dict(zip(hla_df["HLA"], hla_df["protin"]))


def calculate_peptide_features(df: pd.DataFrame, peptide_col: str) -> pd.DataFrame:
    pandas2ri.activate()
    importr("Peptides")
    ro.r(R_FEATURE_CODE)
    calculate_features = ro.globalenv["calculate_features"]

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(df)
        result_r = calculate_features(r_df, peptide_col)
        result_df = ro.conversion.rpy2py(result_r)

    return pd.DataFrame(result_df)


def _calc_features_worker(args):
    chunk, peptide_col = args
    return calculate_peptide_features(chunk, peptide_col)


def _split_dataframe(df: pd.DataFrame, n_chunks: int) -> List[pd.DataFrame]:
    n = len(df)
    if n == 0:
        return []
    n_chunks = max(1, min(int(n_chunks), n))
    base = n // n_chunks
    rem = n % n_chunks
    chunks: List[pd.DataFrame] = []
    start = 0
    for i in range(n_chunks):
        end = start + base + (1 if i < rem else 0)
        chunk = df.iloc[start:end]
        if len(chunk) > 0:
            chunks.append(chunk.copy())
        start = end
    return chunks


def parallel_calculate_features(df: pd.DataFrame, peptide_col: str, num_processes: int) -> pd.DataFrame:
    if len(df) == 0:
        return df.copy()

    n = max(1, min(int(num_processes), len(df)))
    if n == 1:
        return calculate_peptide_features(df, peptide_col)

    chunks = _split_dataframe(df, n)
    with NoDaemonPool(processes=len(chunks)) as pool:
        iterator = pool.imap(_calc_features_worker, [(c, peptide_col) for c in chunks])
        results = list(tqdm(iterator, total=len(chunks), desc="Feature chunks", leave=True))

    return pd.concat(results, ignore_index=True)


def build_feature_table(input_df: pd.DataFrame, peptide_col: str, num_processes: int) -> pd.DataFrame:
    unique_peptides = input_df[[peptide_col]].drop_duplicates().reset_index(drop=True)
    feat = parallel_calculate_features(unique_peptides, peptide_col, num_processes)
    feat = feat.rename(
        columns={
            "aa_composition": "AA Composition",
            "polarity": "Polarity",
            "volume": "Volume",
            "net_charge": "Net Charge",
            "hydrophobicity": "Hydrophobicity",
            "boman_index": "Boman Index",
            "aliphatic_index": "Aliphatic Index",
            "isoelectric_point": "Isoelectric Point",
        }
    )
    return input_df.merge(feat, on=peptide_col, how="left")


def aaindex_vocab(seq: str, vocab_map: Dict[str, int]) -> np.ndarray:
    encoded = np.empty([len(seq)], dtype=np.float32)
    for i, aa in enumerate(seq):
        query = "-" if aa.upper() == "X" else aa.upper()
        encoded[i] = vocab_map.get(query, 0)
    return encoded


def extract_feature_vector(df: pd.DataFrame) -> np.ndarray:
    all_features: List[List[float]] = []
    for i in range(len(df)):
        aa_comp = [float(x.strip()) for x in str(df.iloc[i]["AA Composition"]).split(",")]
        other_feats = df.iloc[i][
            ["Polarity", "Volume", "Net Charge", "Hydrophobicity", "Boman Index", "Aliphatic Index", "Isoelectric Point"]
        ].tolist()
        all_features.append(aa_comp + [float(x) for x in other_feats])
    return np.array(all_features, dtype=np.float32)


def _encode_chunk_worker(args):
    features_chunk, peptides_chunk, hlas_chunk, hla2prot, vocab_map, padding_length, idx_chunk = args
    encoded_batch = []
    flags = []
    out_idx = []

    for i in range(len(peptides_chunk)):
        peptide = peptides_chunk[i]
        hla = hlas_chunk[i]
        hla_paratope = hla2prot.get(hla)

        if hla_paratope is None:
            encoded_batch.append(np.zeros(padding_length, dtype=np.float32))
            flags.append(0)
            out_idx.append(idx_chunk[i])
            continue

        peptide_encode = aaindex_vocab(peptide, vocab_map)
        hla_encode = aaindex_vocab(hla_paratope, vocab_map)
        feature_array = np.asarray(features_chunk[i], dtype=np.float32)

        zero_array = np.zeros(12, dtype=np.float32)
        merged = np.concatenate((feature_array, peptide_encode, zero_array, hla_encode), axis=0)

        if merged.shape[0] > padding_length:
            merged = merged[:padding_length]

        pad_len = padding_length - merged.shape[0]
        merged_padded = np.pad(merged, (0, pad_len), "constant")
        encoded_batch.append(merged_padded.astype(np.float32))
        flags.append(1)
        out_idx.append(idx_chunk[i])

    return encoded_batch, flags, out_idx


def parallel_encode(
    features_np: np.ndarray,
    peptides: List[str],
    hlas: List[str],
    hla2prot: Dict[str, str],
    vocab_map: Dict[str, int],
    padding_length: int,
    num_processes: int,
) -> Tuple[torch.Tensor, List[int]]:
    n = len(peptides)
    if n == 0:
        return torch.empty((0, padding_length), dtype=torch.float32), []

    p = max(1, min(int(num_processes), n))
    if p == 1:
        enc, flags, _ = _encode_chunk_worker(
            (features_np, peptides, hlas, hla2prot, vocab_map, padding_length, list(range(n)))
        )
        return torch.from_numpy(np.stack(enc, axis=0)), flags

    chunk_size = n // p + int(n % p > 0)
    args = []
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        args.append(
            (
                features_np[start:end],
                peptides[start:end],
                hlas[start:end],
                hla2prot,
                vocab_map,
                padding_length,
                list(range(start, end)),
            )
        )

    with NoDaemonPool(processes=len(args)) as pool:
        iterator = pool.imap(_encode_chunk_worker, args)
        results = list(tqdm(iterator, total=len(args), desc="Encoding chunks", leave=True))

    encoded_arrays = [arr for (enc, _, _) in results for arr in enc]
    flags = [f for (_, flg, _) in results for f in flg]
    row_idx = [ix for (_, _, idxs) in results for ix in idxs]

    if row_idx != list(range(n)):
        reorder = np.argsort(np.array(row_idx))
        encoded_arrays = [encoded_arrays[i] for i in reorder]
        flags = [flags[i] for i in reorder]

    encoded_tensor = torch.from_numpy(np.stack(encoded_arrays, axis=0))
    return encoded_tensor, flags


def run_inference(input_df: pd.DataFrame, cfg: InferenceConfig) -> pd.DataFrame:
    if cfg.peptide_col not in input_df.columns or cfg.hla_col not in input_df.columns:
        raise KeyError(f"Input must contain columns '{cfg.peptide_col}' and '{cfg.hla_col}'.")

    t0 = time.time()
    if cfg.verbose:
        print("[Immunogenicity] Start inference")
        print(
            f"[Immunogenicity] rows={len(input_df)} "
            f"batch_size={cfg.batch_size} num_processes={cfg.num_processes} device={cfg.device}"
        )

    work_df = input_df.copy()
    work_df[cfg.peptide_col] = work_df[cfg.peptide_col].astype(str).str.strip().str.upper()
    work_df[cfg.hla_col] = work_df[cfg.hla_col].astype(str).str.strip()
    work_df["_norm_hla"] = work_df[cfg.hla_col].apply(normalize_hla)

    t_feat = time.time()
    if cfg.verbose:
        print("[Immunogenicity] Extracting peptide physicochemical features...")
    work_df = build_feature_table(work_df, cfg.peptide_col, cfg.num_processes)
    if cfg.verbose:
        print(f"[Immunogenicity] Feature extraction done in {time.time() - t_feat:.2f}s")

    infer_df = work_df[
        [
            cfg.peptide_col,
            cfg.hla_col,
            "AA Composition",
            "Polarity",
            "Volume",
            "Net Charge",
            "Hydrophobicity",
            "Boman Index",
            "Aliphatic Index",
            "Isoelectric Point",
        ]
    ].copy()

    infer_df["_norm_hla"] = infer_df[cfg.hla_col].apply(normalize_hla)
    infer_df = infer_df.drop_duplicates(subset=["_norm_hla", cfg.peptide_col]).reset_index(drop=True)

    device = pick_device(cfg.device)
    vocab = Vocab(list(AMINO), min_freq=1)
    hla2prot = parse_hla_fasta(cfg.hla_fasta_path)

    features_np = extract_feature_vector(infer_df)
    peptides = infer_df[cfg.peptide_col].tolist()
    hlas = infer_df["_norm_hla"].tolist()

    t_enc = time.time()
    if cfg.verbose:
        print("[Immunogenicity] Encoding peptide/HLA pairs...")
    encoded_all, flags_all = parallel_encode(
        features_np=features_np,
        peptides=peptides,
        hlas=hlas,
        hla2prot=hla2prot,
        vocab_map=vocab.token_to_idx,
        padding_length=PADDING_LENGTH,
        num_processes=cfg.num_processes,
    )
    if cfg.verbose:
        print(f"[Immunogenicity] Encoding done in {time.time() - t_enc:.2f}s")

    model = build_model(vocab_size=len(vocab), pad_idx=vocab["<pad>"])
    load_model_weights(model, cfg.model_path, device)
    model.to(device)
    model.eval()

    probs_out: List[float] = []
    t_inf = time.time()
    if cfg.verbose:
        print("[Immunogenicity] Running model inference...")
    with torch.no_grad():
        for start in tqdm(range(0, len(infer_df), cfg.batch_size), desc="Infer batches", leave=True):
            end = min(start + cfg.batch_size, len(infer_df))
            batch_encoded = encoded_all[start:end].to(device)
            batch_flags = flags_all[start:end]

            x1 = batch_encoded[:, X2_DIM:].long()
            x2 = batch_encoded[:, 0:X2_DIM].float()
            logits = model(x1, x2)
            probs = torch.softmax(logits, dim=1).cpu().numpy()

            for j, flag in enumerate(batch_flags):
                probs_out.append(float(probs[j, 1]) if flag else 0.5)

    infer_df[cfg.output_score_col] = probs_out

    scored = work_df.merge(
        infer_df[[cfg.peptide_col, "_norm_hla", cfg.output_score_col]],
        on=[cfg.peptide_col, "_norm_hla"],
        how="left",
    )
    scored = scored.drop(columns=["_norm_hla"])

    if cfg.verbose:
        print(f"[Immunogenicity] Inference done in {time.time() - t_inf:.2f}s")
        print(f"[Immunogenicity] Total time {time.time() - t0:.2f}s")

    return scored


def export_model_to_onnx(
    model_path: str,
    onnx_path: str,
    device_name: str = "cpu",
    opset_version: int = 17,
) -> None:
    device = pick_device(device_name)
    vocab = Vocab(list(AMINO), min_freq=1)
    model = build_model(vocab_size=len(vocab), pad_idx=vocab["<pad>"])
    load_model_weights(model, model_path, device)
    model.to(device)
    model.eval()

    x1_len = PADDING_LENGTH - X2_DIM
    dummy_x1 = torch.zeros((1, x1_len), dtype=torch.long, device=device)
    dummy_x2 = torch.zeros((1, X2_DIM), dtype=torch.float32, device=device)

    torch.onnx.export(
        model,
        (dummy_x1, dummy_x2),
        onnx_path,
        export_params=True,
        opset_version=opset_version,
        do_constant_folding=True,
        input_names=["x1", "x2"],
        output_names=["logits"],
        dynamic_axes={
            "x1": {0: "batch", 1: "seq_len"},
            "x2": {0: "batch"},
            "logits": {0: "batch"},
        },
    )
