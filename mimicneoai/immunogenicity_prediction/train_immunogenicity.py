#!/usr/bin/env python3
"""Train MimicNeoAI immunogenicity models from CSV inputs.

This script keeps the original notebook data representation and model path:
25 peptide physicochemical features + peptide AAIndex/PCA indices + 12 zeros
+ HLA AAIndex/PCA indices. By default this reproduces the original full-length
HLA input padded to length 419. NetMHC pseudo-sequence HLA inputs can be used
with --hla-source netmhc-pseudoseq.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
import random
import re
import subprocess
import sys
import time
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score
from torch import nn
from torch.utils import data

try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(iterable=None, *args, **kwargs):
        return iterable if iterable is not None else []


AMINO = "ARNDCQEGHILKMFPSTWYV-"
PADDING_LENGTH = 419
X2_DIM = 25
SEPARATOR_LENGTH = 12
FEATURE_COLUMNS = [
    "aa_composition",
    "polarity",
    "volume",
    "net_charge",
    "hydrophobicity",
    "boman_index",
    "aliphatic_index",
    "isoelectric_point",
]
ENCODING_CACHE_VERSION = "peptide_hla_tensor_v1"


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = SCRIPT_DIR.parent
DEFAULT_MODEL_FILE = SCRIPT_DIR / "model.py"
DEFAULT_FEATURE_RSCRIPT = SCRIPT_DIR / "compute_peptide_features.R"
DEFAULT_AFTER_PCA = SCRIPT_DIR / "resources" / "after_pca.txt"
DEFAULT_HLA_FASTA = PACKAGE_ROOT / "example" / "immunogenicity_prediction" / "models" / "hla_prot.fasta"
DEFAULT_HLA_PSEUDOSEQ_DIR = SCRIPT_DIR / "resources" / "hla_pseudoseq" / "local"
DEFAULT_HLA_CLASS1_PSEUDOSEQ = DEFAULT_HLA_PSEUDOSEQ_DIR / "netmhcpan_class1_allele_to_pseudoseq.csv"
DEFAULT_HLA_CLASS2_PSEUDOSEQ = DEFAULT_HLA_PSEUDOSEQ_DIR / "netmhciipan_class2_allele_to_pseudoseq.csv"


def log_step(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


class Vocab:
    def __init__(self, tokens=None, min_freq=0, reserved_tokens=None):
        if tokens is None:
            tokens = []
        if reserved_tokens is None:
            reserved_tokens = []
        counter = collections.Counter()
        for token in tokens:
            if isinstance(token, list):
                counter.update(token)
            else:
                counter[token] += 1
        self.token_freqs = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        self.idx_to_token = ["<pad>"] + reserved_tokens
        self.token_to_idx = {token: idx for idx, token in enumerate(self.idx_to_token)}
        for token, freq in self.token_freqs:
            if freq >= min_freq and token not in self.token_to_idx:
                self.idx_to_token.append(token)
                self.token_to_idx[token] = len(self.idx_to_token) - 1

    def __len__(self):
        return len(self.idx_to_token)

    def __getitem__(self, tokens):
        if not isinstance(tokens, (list, tuple)):
            return self.token_to_idx.get(tokens, 0)
        return [self.__getitem__(token) for token in tokens]


def load_new_bilstm(model_file: Path):
    spec = importlib.util.spec_from_file_location("mimicneoai_model", model_file)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import model file: {model_file}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NewBiLSTM


def read_table_auto(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


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


def parse_hla_fasta(hla_fasta_path: Path) -> Dict[str, str]:
    hla_names: List[str] = []
    prots: List[str] = []
    prot = ""
    started = False
    with hla_fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
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


def parse_hla_pseudoseq_csv(paths: Sequence[Path]) -> Dict[str, str]:
    hla2seq: Dict[str, str] = {}
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing HLA pseudo-sequence CSV(s): "
            + ", ".join(missing)
            + ". Place local NetMHC-derived resources under "
            + str(DEFAULT_HLA_PSEUDOSEQ_DIR)
            + " or pass --hla-pseudoseq-csv."
        )
    for path in paths:
        df = pd.read_csv(path)
        required = {"allele", "pseudo_sequence"}
        if not required.issubset(df.columns):
            raise KeyError(f"{path} must contain columns: {sorted(required)}")
        for row in df.itertuples(index=False):
            allele = str(getattr(row, "allele")).strip()
            seq = str(getattr(row, "pseudo_sequence")).strip().upper()
            if not allele or not seq or seq.lower() == "nan":
                continue
            keys = {
                allele,
                allele.replace("HLA-", ""),
                normalize_hla(allele),
            }
            for key in keys:
                if key:
                    hla2seq[key] = seq
    return hla2seq


def load_hla_sequences(
    source: str,
    hla_fasta: Path,
    hla_pseudoseq_csv: Sequence[Path] | None,
) -> Dict[str, str]:
    if source == "fasta":
        return parse_hla_fasta(hla_fasta)
    if source == "netmhc-pseudoseq":
        paths = list(hla_pseudoseq_csv or [DEFAULT_HLA_CLASS1_PSEUDOSEQ, DEFAULT_HLA_CLASS2_PSEUDOSEQ])
        return parse_hla_pseudoseq_csv(paths)
    raise ValueError(f"Unsupported HLA source: {source}")


def standardize_input(df: pd.DataFrame, source: str) -> pd.DataFrame:
    rename = {}
    for col in df.columns:
        low = col.strip().lower()
        if low == "peptide":
            rename[col] = "peptide"
        elif low == "hla":
            rename[col] = "hla"
        elif low == "label":
            rename[col] = "label"
        elif low == "length":
            rename[col] = "length"
    df = df.rename(columns=rename).copy()
    missing = {"peptide", "hla", "label"} - set(df.columns)
    if missing:
        raise KeyError(f"{source} missing required columns: {sorted(missing)}")
    df["peptide"] = df["peptide"].astype(str).str.strip().str.upper()
    df["hla"] = df["hla"].astype(str).str.strip()
    df["_norm_hla"] = df["hla"].map(normalize_hla)
    df["label"] = df["label"].astype(int)
    return df


def sha_for_peptides(peptides: Sequence[str]) -> str:
    h = hashlib.sha256()
    for pep in sorted(set(peptides)):
        h.update(pep.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()[:16]


def sha_for_hla_sequences(hla2seq: Dict[str, str]) -> str:
    h = hashlib.sha256()
    for allele, seq in sorted(hla2seq.items()):
        h.update(str(allele).encode("utf-8"))
        h.update(b"\t")
        h.update(str(seq).encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()[:16]


def sha_for_dataframe_columns(df: pd.DataFrame, columns: Sequence[str]) -> str:
    h = hashlib.sha256()
    subset = df.loc[:, list(columns)].copy()
    h.update("\t".join(columns).encode("utf-8"))
    h.update(b"\n")
    h.update(str(len(subset)).encode("utf-8"))
    h.update(b"\n")
    row_hashes = pd.util.hash_pandas_object(subset, index=False).to_numpy(dtype=np.uint64)
    h.update(row_hashes.tobytes())
    return h.hexdigest()[:16]


def ensure_peptide_features(
    df: pd.DataFrame,
    cache_dir: Path,
    rscript_path: Path,
    peptide_col: str = "peptide",
    feature_workers: int = 1,
    feature_chunk_size: int = 250000,
) -> pd.DataFrame:
    cache_dir.mkdir(parents=True, exist_ok=True)
    unique = pd.DataFrame({peptide_col: sorted(df[peptide_col].dropna().astype(str).unique())})
    cache_key = sha_for_peptides(unique[peptide_col].tolist())
    out_csv = cache_dir / f"peptide_features_{cache_key}.csv"
    in_csv = cache_dir / f"peptide_features_{cache_key}.input.csv"
    if not out_csv.exists():
        unique.to_csv(in_csv, index=False)
        t0 = time.time()
        if feature_workers <= 1 or len(unique) <= feature_chunk_size:
            cmd = ["Rscript", str(rscript_path), str(in_csv), str(out_csv), peptide_col]
            log_step(f"[features] running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        else:
            feature_workers = max(1, int(feature_workers))
            feature_chunk_size = max(1, int(feature_chunk_size))
            n_chunks = math.ceil(len(unique) / feature_chunk_size)
            chunks_dir = cache_dir / f"peptide_features_{cache_key}.chunks"
            chunks_dir.mkdir(parents=True, exist_ok=True)
            jobs = []
            for chunk_idx in range(n_chunks):
                start = chunk_idx * feature_chunk_size
                end = min(start + feature_chunk_size, len(unique))
                chunk = unique.iloc[start:end].copy()
                chunk_in = chunks_dir / f"chunk_{chunk_idx:04d}.input.csv"
                chunk_out = chunks_dir / f"chunk_{chunk_idx:04d}.features.csv"
                chunk.to_csv(chunk_in, index=False)
                jobs.append((chunk_idx, chunk_in, chunk_out, len(chunk)))

            log_step(
                f"[features] running parallel Rscript workers={feature_workers} "
                f"chunks={n_chunks} chunk_size={feature_chunk_size} rows={len(unique)}"
            )
            active: list[tuple[int, subprocess.Popen, Path, int]] = []
            completed = 0
            next_job = 0
            while next_job < len(jobs) or active:
                while next_job < len(jobs) and len(active) < feature_workers:
                    chunk_idx, chunk_in, chunk_out, n_rows = jobs[next_job]
                    cmd = ["Rscript", str(rscript_path), str(chunk_in), str(chunk_out), peptide_col]
                    log_step(f"[features] start chunk={chunk_idx + 1}/{n_chunks} rows={n_rows}")
                    proc = subprocess.Popen(cmd)
                    active.append((chunk_idx, proc, chunk_out, n_rows))
                    next_job += 1

                still_active = []
                for chunk_idx, proc, chunk_out, n_rows in active:
                    rc = proc.poll()
                    if rc is None:
                        still_active.append((chunk_idx, proc, chunk_out, n_rows))
                        continue
                    if rc != 0:
                        for _, other_proc, _, _ in still_active:
                            other_proc.terminate()
                        raise subprocess.CalledProcessError(rc, proc.args)
                    if not chunk_out.exists():
                        raise FileNotFoundError(f"Missing peptide feature chunk output: {chunk_out}")
                    completed += 1
                    log_step(f"[features] done chunk={chunk_idx + 1}/{n_chunks} completed={completed}/{n_chunks}")
                active = still_active
                if active:
                    time.sleep(1)

            header = True
            with out_csv.open("w", encoding="utf-8", newline="") as out_handle:
                for chunk_idx, _, chunk_out, _ in jobs:
                    with chunk_out.open("r", encoding="utf-8") as in_handle:
                        for line_no, line in enumerate(in_handle):
                            if line_no == 0:
                                if header:
                                    out_handle.write(line)
                                    header = False
                                continue
                            out_handle.write(line)
        log_step(f"[features] done: rows={len(unique)} elapsed={time.time() - t0:.1f}s")
    else:
        log_step(f"[features] cache hit: {out_csv.name} rows={len(unique)}")
    feat = pd.read_csv(out_csv)
    return df.merge(feat[[peptide_col] + FEATURE_COLUMNS], on=peptide_col, how="left")


def aaindex_vocab(seq: str, vocab_map: Dict[str, int]) -> np.ndarray:
    encoded = np.empty([len(seq)], dtype=np.float32)
    for i, aa in enumerate(seq):
        query = "-" if aa.upper() == "X" else aa.upper()
        encoded[i] = vocab_map.get(query, 0)
    return encoded


def extract_feature_vector(row: pd.Series) -> List[float]:
    aa_comp = [float(value.strip()) for value in str(row["aa_composition"]).split(",")]
    other = [float(row[col]) for col in FEATURE_COLUMNS[1:]]
    features = aa_comp + other
    if len(features) != X2_DIM:
        raise ValueError(f"Expected {X2_DIM} peptide features, got {len(features)}")
    return features


def resolve_padding_length(
    requested_padding_length: int,
    hla_source: str,
    combined_df: pd.DataFrame,
    hla2seq: Dict[str, str],
) -> int:
    if requested_padding_length > 0:
        return int(requested_padding_length)
    if hla_source == "fasta":
        return PADDING_LENGTH
    max_peptide_len = int(combined_df["peptide"].astype(str).str.len().max())
    max_hla_len = max(len(seq) for seq in hla2seq.values())
    return X2_DIM + max_peptide_len + SEPARATOR_LENGTH + max_hla_len


def _encode_dataframe_chunk(args):
    chunk_df, hla2prot, vocab_map, padding_length = args
    encoded_rows: List[np.ndarray] = []
    labels: List[int] = []
    kept_indices: List[int] = []
    missing_hla = 0

    for idx, row in chunk_df.iterrows():
        hla_paratope = hla2prot.get(str(row["_norm_hla"]))
        if hla_paratope is None:
            missing_hla += 1
            continue
        feature_array = np.asarray(extract_feature_vector(row), dtype=np.float32)
        peptide_encode = aaindex_vocab(str(row["peptide"]), vocab_map)
        hla_encode = aaindex_vocab(hla_paratope, vocab_map)
        zero_array = np.zeros(SEPARATOR_LENGTH, dtype=np.float32)
        merged = np.concatenate((feature_array, peptide_encode, zero_array, hla_encode), axis=0)
        if merged.shape[0] > padding_length:
            merged = merged[:padding_length]
        pad_len = padding_length - merged.shape[0]
        encoded_rows.append(np.pad(merged, (0, pad_len), "constant").astype(np.float32))
        labels.append(int(row["label"]))
        kept_indices.append(int(idx))

    if encoded_rows:
        encoded = np.stack(encoded_rows, axis=0)
    else:
        encoded = np.empty((0, padding_length), dtype=np.float32)
    return encoded, labels, kept_indices, missing_hla


def _split_dataframe_rows(df: pd.DataFrame, chunk_size: int) -> List[pd.DataFrame]:
    if len(df) == 0:
        return []
    chunk_size = max(1, int(chunk_size))
    return [df.iloc[start : start + chunk_size].copy() for start in range(0, len(df), chunk_size)]


def encode_dataset(
    df: pd.DataFrame,
    hla2prot: Dict[str, str],
    vocab: Vocab,
    padding_length: int,
    encode_workers: int = 1,
    encode_chunk_size: int = 100000,
) -> Tuple[torch.Tensor, torch.Tensor, pd.DataFrame]:
    reset_df = df.reset_index(drop=True)
    chunks = _split_dataframe_rows(reset_df, encode_chunk_size)
    log_step(
        f"[encode] rows={len(reset_df)} chunks={len(chunks)} "
        f"workers={encode_workers} chunk_size={encode_chunk_size}"
    )
    worker_args = [(chunk, hla2prot, vocab.token_to_idx, padding_length) for chunk in chunks]

    encoded_parts: List[np.ndarray] = []
    labels: List[int] = []
    kept_indices: List[int] = []
    missing_hla = 0
    if encode_workers <= 1 or len(chunks) <= 1:
        iterator = tqdm(worker_args, total=len(worker_args), desc="Encoding peptide-HLA chunks", leave=False)
        for args in iterator:
            encoded, chunk_labels, chunk_indices, chunk_missing = _encode_dataframe_chunk(args)
            if len(encoded):
                encoded_parts.append(encoded)
            labels.extend(chunk_labels)
            kept_indices.extend(chunk_indices)
            missing_hla += chunk_missing
    else:
        max_workers = max(1, int(encode_workers))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(_encode_dataframe_chunk, args) for args in worker_args]
            for future in tqdm(
                futures,
                total=len(futures),
                desc=f"Encoding peptide-HLA chunks ({max_workers} workers)",
                leave=False,
            ):
                encoded, chunk_labels, chunk_indices, chunk_missing = future.result()
                if len(encoded):
                    encoded_parts.append(encoded)
                labels.extend(chunk_labels)
                kept_indices.extend(chunk_indices)
                missing_hla += chunk_missing

    if missing_hla:
        log_step(f"[encode] skipped rows with missing HLA sequence: {missing_hla}")
    if not encoded_parts:
        raise RuntimeError("No rows could be encoded. Check HLA fasta coverage.")
    encoded = torch.from_numpy(np.concatenate(encoded_parts, axis=0))
    y = torch.tensor(labels, dtype=torch.long)
    kept = df.reset_index(drop=True).iloc[kept_indices].reset_index(drop=True)
    kept.attrs["_kept_indices"] = kept_indices
    return encoded, y, kept


def encode_dataset_cached(
    df: pd.DataFrame,
    hla2prot: Dict[str, str],
    vocab: Vocab,
    padding_length: int,
    cache_dir: Path,
    split_name: str,
    hla_source: str,
    hla_digest: str,
    enabled: bool = True,
    encode_workers: int = 1,
    encode_chunk_size: int = 100000,
) -> Tuple[torch.Tensor, torch.Tensor, pd.DataFrame]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    encoded_columns = ["peptide", "_norm_hla", "label"] + FEATURE_COLUMNS
    data_digest = sha_for_dataframe_columns(df, encoded_columns)
    cache_key_payload = {
        "version": ENCODING_CACHE_VERSION,
        "split": split_name,
        "data_digest": data_digest,
        "hla_source": hla_source,
        "hla_digest": hla_digest,
        "padding_length": int(padding_length),
        "vocab": AMINO,
        "x2_dim": X2_DIM,
        "separator_length": SEPARATOR_LENGTH,
    }
    key_hash = hashlib.sha256(
        json.dumps(cache_key_payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    cache_path = cache_dir / f"encoded_{split_name}_{key_hash}.pt"

    reset_df = df.reset_index(drop=True)
    if enabled and cache_path.exists():
        try:
            payload = torch.load(cache_path, map_location="cpu", weights_only=False)
        except TypeError:
            payload = torch.load(cache_path, map_location="cpu")
        metadata = payload.get("metadata", {})
        if metadata == cache_key_payload:
            kept_indices = payload["kept_indices"]
            kept = reset_df.iloc[kept_indices].reset_index(drop=True)
            log_step(
                f"[encode] cache hit split={split_name} file={cache_path.name} rows={len(kept_indices)}"
            )
            return payload["x"], payload["y"], kept
        log_step(f"[encode] cache metadata mismatch, rebuilding split={split_name}: {cache_path.name}")

    x, y, kept = encode_dataset(
        reset_df,
        hla2prot,
        vocab,
        padding_length,
        encode_workers=encode_workers,
        encode_chunk_size=encode_chunk_size,
    )
    if enabled:
        kept_indices = [int(i) for i in kept.attrs.get("_kept_indices", kept.index.tolist())]
        torch.save(
            {
                "x": x.cpu(),
                "y": y.cpu(),
                "kept_indices": kept_indices,
                "metadata": cache_key_payload,
            },
            cache_path,
        )
        log_step(f"[encode] cache write split={split_name} file={cache_path.name} rows={len(kept_indices)}")
    return x, y, kept


def initialize_model(new_bilstm, after_pca_path: Path, device: torch.device) -> nn.Module:
    after_pca = np.loadtxt(after_pca_path)
    zero_row = np.zeros((1, after_pca.shape[1]), dtype=np.float32)
    embeds = torch.tensor(np.vstack([zero_row, after_pca]), dtype=torch.float32)
    vocab = Vocab(list(AMINO), min_freq=1)
    model = new_bilstm(
        vocab_size=len(vocab),
        embedding_dim=12,
        hidden_dim_x1=256,
        hidden_dim_x2=16,
        output_dim=2,
        n_layers=2,
        bidirectional=True,
        dropout=0.5,
        pad_idx=vocab["<pad>"],
        x2_dim=25,
    )
    model.embedding.weight.data.copy_(embeds)
    model.embedding.weight.requires_grad = False
    return model.to(device)


def load_initial_weights(model: nn.Module, model_path: Path, device: torch.device) -> None:
    state = torch.load(model_path, map_location=device)
    if isinstance(state, dict) and "model_state" in state:
        state = state["model_state"]
    if not isinstance(state, dict):
        raise RuntimeError(f"Unsupported checkpoint format: {model_path}")

    candidates = [state]

    stripped = OrderedDict()
    for key, value in state.items():
        stripped[key.replace("module.", "", 1)] = value
    candidates.append(stripped)

    prefixed = OrderedDict()
    for key, value in state.items():
        prefixed[key if key.startswith("module.") else f"module.{key}"] = value
    candidates.append(prefixed)

    last_error = None
    for candidate in candidates:
        try:
            model.load_state_dict(candidate, strict=True)
            return
        except RuntimeError as exc:
            last_error = exc
    raise RuntimeError(f"Failed to load initial weights from {model_path}. Last error: {last_error}")


def accuracy(logits: torch.Tensor, labels: torch.Tensor) -> float:
    return float((logits.argmax(dim=1) == labels).sum().item())


def evaluate(model: nn.Module, loader: data.DataLoader, device: torch.device) -> Dict[str, float]:
    model.eval()
    labels_all: List[int] = []
    probs_all: List[float] = []
    pred_all: List[int] = []
    correct = 0.0
    total = 0
    with torch.no_grad():
        for features, labels in tqdm(loader, desc="Evaluating", leave=False):
            features = features.to(device)
            labels = labels.to(device)
            x1 = features[:, X2_DIM:].long()
            x2 = features[:, 0:X2_DIM].float()
            logits = model(x1, x2)
            probs = torch.softmax(logits, dim=1)[:, 1]
            correct += accuracy(logits, labels)
            total += labels.numel()
            labels_all.extend(labels.detach().cpu().numpy().astype(int).tolist())
            probs_all.extend(probs.detach().cpu().numpy().astype(float).tolist())
            pred_all.extend(logits.argmax(dim=1).detach().cpu().numpy().astype(int).tolist())
    metrics = {"accuracy": correct / max(total, 1)}
    if len(set(labels_all)) == 2:
        metrics["roc_auc"] = float(roc_auc_score(labels_all, probs_all))
        metrics["average_precision"] = float(average_precision_score(labels_all, probs_all))
    else:
        metrics["roc_auc"] = float("nan")
        metrics["average_precision"] = float("nan")
    return metrics


def predict(model: nn.Module, loader: data.DataLoader, device: torch.device) -> Tuple[List[int], List[float], List[int]]:
    model.eval()
    labels_all: List[int] = []
    probs_all: List[float] = []
    pred_all: List[int] = []
    with torch.no_grad():
        for features, labels in tqdm(loader, desc="Predicting", leave=False):
            features = features.to(device)
            labels = labels.to(device)
            x1 = features[:, X2_DIM:].long()
            x2 = features[:, 0:X2_DIM].float()
            logits = model(x1, x2)
            probs = torch.softmax(logits, dim=1)[:, 1]
            labels_all.extend(labels.detach().cpu().numpy().astype(int).tolist())
            probs_all.extend(probs.detach().cpu().numpy().astype(float).tolist())
            pred_all.extend(logits.argmax(dim=1).detach().cpu().numpy().astype(int).tolist())
    return labels_all, probs_all, pred_all


def train_one_model(
    model: nn.Module,
    train_loader: data.DataLoader,
    val_loader: data.DataLoader | None,
    device: torch.device,
    epochs: int,
    lr: float,
    checkpoint_dir: Path,
    model_name: str,
    loss_reduction: str,
    class_weights: torch.Tensor | None = None,
) -> Tuple[nn.Module, List[Dict[str, float]], Dict[str, Path]]:
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    trainer = torch.optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.CrossEntropyLoss(
        reduction="none",
        weight=class_weights.to(device) if class_weights is not None else None,
    )
    history: List[Dict[str, float]] = []
    best_values = {
        "val_roc_auc": -float("inf"),
        "val_average_precision": -float("inf"),
    }
    best_paths = {
        "val_roc_auc": checkpoint_dir / f"{model_name}_best_val_roc_auc.pth",
        "val_average_precision": checkpoint_dir / f"{model_name}_best_val_average_precision.pth",
    }
    for epoch in range(epochs):
        model.train()
        loss_sum = 0.0
        correct = 0.0
        n = 0
        batch_iter = tqdm(
            train_loader,
            desc=f"{model_name} epoch {epoch + 1}/{epochs}",
            leave=False,
        )
        for features, labels in batch_iter:
            features = features.to(device)
            labels = labels.to(device)
            x1 = features[:, X2_DIM:].long()
            x2 = features[:, 0:X2_DIM].float()
            trainer.zero_grad()
            logits = model(x1, x2)
            loss = loss_fn(logits, labels)
            if loss_reduction == "mean":
                loss.mean().backward()
            else:
                loss.sum().backward()
            trainer.step()
            loss_sum += float(loss.sum().detach().cpu().item())
            correct += accuracy(logits, labels)
            n += labels.numel()
            batch_iter.set_postfix(
                loss=f"{loss_sum / max(n, 1):.4f}",
                acc=f"{correct / max(n, 1):.4f}",
            )
        row = {
            "epoch": epoch,
            "train_loss": loss_sum / max(n, 1),
            "train_accuracy": correct / max(n, 1),
        }
        if val_loader is not None:
            row.update({f"val_{k}": v for k, v in evaluate(model, val_loader, device).items()})
            for metric_name, best_path in best_paths.items():
                metric_value = float(row.get(metric_name, float("nan")))
                if np.isfinite(metric_value) and metric_value > best_values[metric_name]:
                    best_values[metric_name] = metric_value
                    torch.save(
                        {
                            "epoch": epoch,
                            "model_state": model.state_dict(),
                            "optimizer_state": trainer.state_dict(),
                            "history": history + [row],
                            "best_metric": metric_name,
                            "best_value": metric_value,
                        },
                        best_path,
                    )
        history.append(row)
        msg = (
            f"[{model_name}] epoch={epoch + 1}/{epochs} "
            f"loss={row['train_loss']:.4f} train_acc={row['train_accuracy']:.4f}"
        )
        if val_loader is not None:
            msg += (
                f" val_acc={row.get('val_accuracy', float('nan')):.4f}"
                f" val_auc={row.get('val_roc_auc', float('nan')):.4f}"
                f" val_ap={row.get('val_average_precision', float('nan')):.4f}"
            )
        log_step(msg)
        if (epoch + 1) % 100 == 0 or epoch == epochs - 1:
            torch.save(
                {
                    "epoch": epoch,
                    "model_state": model.state_dict(),
                    "optimizer_state": trainer.state_dict(),
                    "history": history,
                },
                checkpoint_dir / f"{model_name}_epoch_{epoch}.pth",
            )
    return model, history, best_paths


def old_notebook_kfold_indices(labels: Sequence[int], k: int) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    labels_np = np.asarray(labels)
    class_indices = [np.where(labels_np == cls)[0].tolist() for cls in sorted(set(labels_np.tolist()))]
    fold_chunks: List[List[List[int]]] = []
    for idxs in class_indices:
        chunk_size = int(len(idxs) / k)
        chunks = []
        start = 0
        for fold in range(k):
            if fold != k - 1:
                chunks.append(idxs[start : start + chunk_size])
                start += chunk_size
            else:
                chunks.append(idxs[start:])
        fold_chunks.append(chunks)
    for fold in range(k):
        test_idx: List[int] = []
        train_idx: List[int] = []
        for chunks in fold_chunks:
            for i, chunk in enumerate(chunks):
                if i == fold:
                    test_idx.extend(chunk)
                else:
                    train_idx.extend(chunk)
        yield np.asarray(train_idx, dtype=int), np.asarray(test_idx, dtype=int)


def make_loader(features: torch.Tensor, labels: torch.Tensor, indices: np.ndarray, batch_size: int, shuffle: bool) -> data.DataLoader:
    ds = data.TensorDataset(features[indices], labels[indices])
    return data.DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=False)


def make_training_loader(
    features: torch.Tensor,
    labels: torch.Tensor,
    kept_df: pd.DataFrame,
    batch_size: int,
    sampler_mode: str,
    sampler_source_col: str,
    outdir: Path,
) -> data.DataLoader:
    ds = data.TensorDataset(features, labels)
    if sampler_mode == "shuffle":
        return data.DataLoader(ds, batch_size=batch_size, shuffle=True, drop_last=False)

    meta = kept_df.reset_index(drop=True).copy()
    if len(meta) != len(labels):
        raise ValueError(f"Sampler metadata rows ({len(meta)}) != labels ({len(labels)})")
    meta["_label"] = labels.detach().cpu().numpy().astype(int)
    if sampler_mode == "label_balanced":
        group_cols = ["_label"]
    elif sampler_mode == "source_label_balanced":
        if sampler_source_col not in meta.columns:
            raise KeyError(f"Missing sampler source column: {sampler_source_col}")
        group_cols = [sampler_source_col, "_label"]
    else:
        raise ValueError(f"Unsupported train sampler: {sampler_mode}")

    group_sizes = meta.groupby(group_cols, dropna=False).size().rename("n").reset_index()
    meta = meta.merge(group_sizes, on=group_cols, how="left")
    weights = (1.0 / meta["n"].astype(float)).to_numpy(dtype=np.float64)
    sampler = data.WeightedRandomSampler(
        weights=torch.as_tensor(weights, dtype=torch.double),
        num_samples=len(weights),
        replacement=True,
    )
    report_path = outdir / "metrics" / "train_sampler_groups.tsv"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    group_sizes.to_csv(report_path, sep="\t", index=False)
    log_step(
        f"[sampler] mode={sampler_mode} groups={len(group_sizes)} "
        f"source_col={sampler_source_col} report={report_path}"
    )
    return data.DataLoader(ds, batch_size=batch_size, sampler=sampler, drop_last=False)


def compute_class_weights(labels: torch.Tensor, mode: str) -> torch.Tensor | None:
    if mode == "none":
        return None
    if mode != "balanced":
        raise ValueError(f"Unsupported class weight mode: {mode}")
    label_np = labels.detach().cpu().numpy().astype(int)
    counts = np.bincount(label_np, minlength=2).astype(np.float64)
    if np.any(counts == 0):
        raise ValueError(f"Cannot compute balanced class weights with empty class: {counts.tolist()}")
    weights = counts.sum() / (len(counts) * counts)
    out = torch.tensor(weights, dtype=torch.float32)
    log_step(f"[class_weight] mode=balanced counts={counts.astype(int).tolist()} weights={weights.tolist()}")
    return out


def write_metrics(path: Path, rows: List[Dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--antigen", required=True)
    parser.add_argument("--train-csv", "--train-table", dest="train_csv", required=True, type=Path)
    parser.add_argument("--val-csv", "--val-table", dest="val_csv", required=True, type=Path)
    parser.add_argument("--test-csv", "--test-table", dest="test_csv", default=None, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--model-file", default=DEFAULT_MODEL_FILE, type=Path)
    parser.add_argument(
        "--hla-source",
        choices=["fasta", "netmhc-pseudoseq"],
        default="fasta",
        help="Use original full-length HLA fasta or local NetMHC pseudo-sequence CSV resources.",
    )
    parser.add_argument("--hla-fasta", default=DEFAULT_HLA_FASTA, type=Path)
    parser.add_argument(
        "--hla-pseudoseq-csv",
        action="append",
        default=None,
        type=Path,
        help="NetMHC pseudo-sequence CSV with allele,pseudo_sequence columns. Can be passed multiple times.",
    )
    parser.add_argument(
        "--padding-length",
        default=0,
        type=int,
        help="Sequence padding length. 0 keeps 419 for fasta and auto-sizes for NetMHC pseudo-sequences.",
    )
    parser.add_argument("--after-pca", default=DEFAULT_AFTER_PCA, type=Path)
    parser.add_argument("--feature-rscript", default=DEFAULT_FEATURE_RSCRIPT, type=Path)
    parser.add_argument("--feature-workers", default=1, type=int)
    parser.add_argument("--feature-chunk-size", default=250000, type=int)
    parser.add_argument("--encode-workers", default=1, type=int)
    parser.add_argument("--encode-chunk-size", default=100000, type=int)
    parser.add_argument(
        "--no-encoded-cache",
        dest="encoded_cache",
        action="store_false",
        help="Disable cached encoded tensor reuse.",
    )
    parser.set_defaults(encoded_cache=True)
    parser.add_argument("--init-from", default=None, type=Path)
    parser.add_argument(
        "--select-best-metric",
        choices=["last", "val_roc_auc", "val_average_precision"],
        default="last",
        help="Checkpoint used for final validation/test prediction after training.",
    )
    parser.add_argument("--epochs", default=200, type=int)
    parser.add_argument("--batch-size", default=256, type=int)
    parser.add_argument("--lr", default=5e-5, type=float)
    parser.add_argument(
        "--loss-reduction",
        choices=["sum", "mean"],
        default="sum",
        help="Use sum to reproduce the original notebooks; use mean for batch-size-stable tuning.",
    )
    parser.add_argument(
        "--class-weight-mode",
        choices=["none", "balanced"],
        default="none",
        help="Optional class weights for CrossEntropyLoss.",
    )
    parser.add_argument(
        "--train-sampler",
        choices=["shuffle", "label_balanced", "source_label_balanced"],
        default="shuffle",
        help="Sampling strategy for the final training loader.",
    )
    parser.add_argument(
        "--sampler-source-col",
        default="dataset_source",
        help="Metadata column used by --train-sampler source_label_balanced.",
    )
    parser.add_argument("--kfold", default=5, type=int)
    parser.add_argument("--seed", default=2026, type=int)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--smoke-rows-per-label", default=0, type=int)
    parser.add_argument("--skip-kfold", action="store_true")
    args = parser.parse_args()

    t0 = time.time()
    set_seed(args.seed)
    log_step(f"[start] antigen={args.antigen} epochs={args.epochs} batch_size={args.batch_size}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for sub in ["cache", "models", "metrics", "predictions"]:
        (args.outdir / sub).mkdir(exist_ok=True)
    config = vars(args).copy()
    config = {k: str(v) if isinstance(v, Path) else v for k, v in config.items()}
    (args.outdir / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    device = torch.device("cuda:0" if args.device == "auto" and torch.cuda.is_available() else args.device)
    log_step(f"[env] device={device} cuda_available={torch.cuda.is_available()}")
    NewBiLSTM = load_new_bilstm(args.model_file)
    vocab = Vocab(list(AMINO), min_freq=1)
    hla2prot = load_hla_sequences(args.hla_source, args.hla_fasta, args.hla_pseudoseq_csv)
    hla_digest = sha_for_hla_sequences(hla2prot)
    hla_lengths = sorted({len(seq) for seq in hla2prot.values()})
    log_step(
        f"[resources] hla_source={args.hla_source} HLA alleles={len(hla2prot)} "
        f"hla_lengths={hla_lengths[:10]} hla_digest={hla_digest} after_pca={args.after_pca}"
    )

    train_df = standardize_input(read_table_auto(args.train_csv), "train")
    val_df = standardize_input(read_table_auto(args.val_csv), "validation")
    test_df = standardize_input(read_table_auto(args.test_csv), "test") if args.test_csv else None
    if args.smoke_rows_per_label > 0:
        train_df = (
            train_df.groupby("label", group_keys=False)
            .head(args.smoke_rows_per_label)
            .reset_index(drop=True)
        )
        val_df = (
            val_df.groupby("label", group_keys=False)
            .head(min(args.smoke_rows_per_label, 10))
            .reset_index(drop=True)
        )
        if test_df is not None:
            test_df = (
                test_df.groupby("label", group_keys=False)
                .head(min(args.smoke_rows_per_label, 10))
                .reset_index(drop=True)
            )

    split_frames = [train_df.assign(_split="train"), val_df.assign(_split="validation")]
    if test_df is not None:
        split_frames.append(test_df.assign(_split="test"))
    combined = pd.concat(split_frames, ignore_index=True)
    combined = ensure_peptide_features(
        combined,
        args.outdir / "cache",
        args.feature_rscript,
        feature_workers=args.feature_workers,
        feature_chunk_size=args.feature_chunk_size,
    )
    padding_length = resolve_padding_length(args.padding_length, args.hla_source, combined, hla2prot)
    log_step(f"[encode] padding_length={padding_length}")
    train_feat = combined[combined["_split"] == "train"].drop(columns=["_split"]).reset_index(drop=True)
    val_feat = combined[combined["_split"] == "validation"].drop(columns=["_split"]).reset_index(drop=True)
    test_feat = (
        combined[combined["_split"] == "test"].drop(columns=["_split"]).reset_index(drop=True)
        if test_df is not None
        else None
    )

    log_step("[encode] training set")
    train_x, train_y, train_kept = encode_dataset_cached(
        train_feat,
        hla2prot,
        vocab,
        padding_length,
        args.outdir / "cache",
        "train",
        args.hla_source,
        hla_digest,
        enabled=args.encoded_cache,
        encode_workers=args.encode_workers,
        encode_chunk_size=args.encode_chunk_size,
    )
    log_step("[encode] validation set")
    val_x, val_y, val_kept = encode_dataset_cached(
        val_feat,
        hla2prot,
        vocab,
        padding_length,
        args.outdir / "cache",
        "validation",
        args.hla_source,
        hla_digest,
        enabled=args.encoded_cache,
        encode_workers=args.encode_workers,
        encode_chunk_size=args.encode_chunk_size,
    )
    if test_feat is not None:
        log_step("[encode] test set")
        test_x, test_y, test_kept = encode_dataset_cached(
            test_feat,
            hla2prot,
            vocab,
            padding_length,
            args.outdir / "cache",
            "test",
            args.hla_source,
            hla_digest,
            enabled=args.encoded_cache,
            encode_workers=args.encode_workers,
            encode_chunk_size=args.encode_chunk_size,
        )
    else:
        test_x, test_y, test_kept = None, None, None
    log_step(
        f"[data] train rows={len(train_df)} encoded={len(train_kept)} "
        f"validation rows={len(val_df)} encoded={len(val_kept)} device={device}"
    )
    class_weights = compute_class_weights(train_y, args.class_weight_mode)

    all_fold_rows: List[Dict[str, float]] = []
    if not args.skip_kfold:
        for fold, (train_idx, test_idx) in enumerate(old_notebook_kfold_indices(train_y.numpy(), args.kfold)):
            log_step(f"[kfold] fold={fold} train={len(train_idx)} test={len(test_idx)}")
            model = initialize_model(NewBiLSTM, args.after_pca, device)
            if args.init_from is not None:
                load_initial_weights(model, args.init_from, device)
            train_loader = make_loader(train_x, train_y, train_idx, args.batch_size, shuffle=True)
            test_loader = make_loader(train_x, train_y, test_idx, args.batch_size, shuffle=False)
            model, history, _ = train_one_model(
                model,
                train_loader,
                test_loader,
                device,
                args.epochs,
                args.lr,
                args.outdir / "models" / f"fold_{fold}",
                f"{args.antigen}_fold_{fold}",
                args.loss_reduction,
                class_weights,
            )
            final_metrics = evaluate(model, test_loader, device)
            final_metrics.update({"fold": fold, "n_train": len(train_idx), "n_test": len(test_idx)})
            all_fold_rows.append(final_metrics)
            write_metrics(args.outdir / "metrics" / f"fold_{fold}_history.tsv", history)
        write_metrics(args.outdir / "metrics" / "kfold_metrics.tsv", all_fold_rows)

    final_model = initialize_model(NewBiLSTM, args.after_pca, device)
    if args.init_from is not None:
        log_step(f"[init] loading weights from {args.init_from}")
        load_initial_weights(final_model, args.init_from, device)
    final_train_loader = make_training_loader(
        train_x,
        train_y,
        train_kept,
        args.batch_size,
        args.train_sampler,
        args.sampler_source_col,
        args.outdir,
    )
    val_loader = data.DataLoader(
        data.TensorDataset(val_x, val_y),
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )
    final_model, final_history, final_best_paths = train_one_model(
        final_model,
        final_train_loader,
        val_loader,
        device,
        args.epochs,
        args.lr,
        args.outdir / "models" / "final_checkpoints",
        f"{args.antigen}_final",
        args.loss_reduction,
        class_weights,
    )
    selected_checkpoint = None
    if args.select_best_metric != "last":
        selected_checkpoint = final_best_paths[args.select_best_metric]
        log_step(f"[select] loading best checkpoint by {args.select_best_metric}: {selected_checkpoint}")
        load_initial_weights(final_model, selected_checkpoint, device)
    torch.save(final_model.state_dict(), args.outdir / "models" / f"MimicNeoAI_{args.antigen}_final.pth")
    write_metrics(args.outdir / "metrics" / "final_history.tsv", final_history)

    labels, probs, preds = predict(final_model, val_loader, device)
    pred_df = val_kept.copy()
    pred_df["true_label"] = labels
    pred_df["pred_label"] = preds
    pred_df["immunogenicity_score"] = probs
    pred_df.to_csv(args.outdir / "predictions" / "validation_predictions.tsv", sep="\t", index=False)
    final_metrics = {
        "accuracy": float(accuracy_score(labels, preds)),
        "roc_auc": float(roc_auc_score(labels, probs)) if len(set(labels)) == 2 else float("nan"),
        "average_precision": float(average_precision_score(labels, probs)) if len(set(labels)) == 2 else float("nan"),
        "n_validation": len(labels),
        "elapsed_seconds": time.time() - t0,
        "selected_checkpoint": str(selected_checkpoint) if selected_checkpoint is not None else "last",
    }
    write_metrics(args.outdir / "metrics" / "validation_metrics.tsv", [final_metrics])
    if test_x is not None and test_y is not None and test_kept is not None:
        test_loader = data.DataLoader(
            data.TensorDataset(test_x, test_y),
            batch_size=args.batch_size,
            shuffle=False,
            drop_last=False,
        )
        test_labels, test_probs, test_preds = predict(final_model, test_loader, device)
        test_pred_df = test_kept.copy()
        test_pred_df["true_label"] = test_labels
        test_pred_df["pred_label"] = test_preds
        test_pred_df["immunogenicity_score"] = test_probs
        test_pred_df.to_csv(args.outdir / "predictions" / "test_predictions.tsv", sep="\t", index=False)
        test_metrics = {
            "accuracy": float(accuracy_score(test_labels, test_preds)),
            "roc_auc": float(roc_auc_score(test_labels, test_probs)) if len(set(test_labels)) == 2 else float("nan"),
            "average_precision": float(average_precision_score(test_labels, test_probs)) if len(set(test_labels)) == 2 else float("nan"),
            "n_test": len(test_labels),
            "elapsed_seconds": time.time() - t0,
            "selected_checkpoint": str(selected_checkpoint) if selected_checkpoint is not None else "last",
        }
        write_metrics(args.outdir / "metrics" / "test_metrics.tsv", [test_metrics])
        log_step("[test] " + json.dumps(test_metrics, indent=2))
    log_step("[done] " + json.dumps(final_metrics, indent=2))


if __name__ == "__main__":
    main()
