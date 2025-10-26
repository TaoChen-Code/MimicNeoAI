# -*- coding: utf-8 -*-

# ===== Imports =====
import collections
import re
import warnings
from collections import OrderedDict
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import nn
from torchvision import models  # kept if used elsewhere
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)  # suppress FutureWarning noise


# ===== Model =====
class newBiLSTM(nn.Module):
    """
    BiLSTM branch for sequence (x1) + FC branch for numeric features (x2); outputs 2-class logits.
    """
    def __init__(
        self,
        vocab_size,
        embedding_dim,
        hidden_dim_x1,
        hidden_dim_x2,
        output_dim,
        n_layers,
        bidirectional,
        dropout,
        pad_idx,
        x2_dim,
    ):
        super(newBiLSTM, self).__init__()

        # Embedding and LSTM for x1 (sequence indices)
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_idx)
        self.lstm = nn.LSTM(
            embedding_dim,
            hidden_dim_x1,
            num_layers=n_layers,
            bidirectional=bidirectional,
            dropout=dropout,
        )

        # FC for x2 (handcrafted features)
        self.fc_x2 = nn.Linear(x2_dim, hidden_dim_x2)

        # Final classifier on concatenated [x1_lstm, x2_fc]
        self.fc = nn.Linear(hidden_dim_x1 * 2 + hidden_dim_x2, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x1, x2):
        # Ensure contiguous params for speed/compatibility on some backends
        self.lstm.flatten_parameters()

        # x1: (batch, seq_len) int indices -> (seq_len, batch, emb)
        embedded = self.dropout(self.embedding(x1.T))
        _, (hidden, _) = self.lstm(embedded)

        # Concatenate last forward/backward hidden states -> (batch, hidden_dim_x1*2)
        hidden_x1 = self.dropout(torch.cat((hidden[-2, :, :], hidden[-1, :, :]), dim=1))

        # x2: (batch, x2_dim)
        hidden_x2 = self.dropout(self.fc_x2(x2))

        # Concatenate and classify
        combined = torch.cat((hidden_x1, hidden_x2), dim=1)
        return self.fc(combined)


# ===== Data helpers =====
def load_hla(hla_prot_file):
    """
    Load HLA -> paratope/protein sequence mapping from a FASTA-like file.

    Header format assumption: first '>' line has a name; we take the second token as HLA name.
    Sequences are concatenated per record. HLA names are normalized to first two fields (e.g., A*02:01).
    """
    hla_names = []
    prots = []
    with open(hla_prot_file, "r") as f:
        lines = f.readlines()
        switch = 0
        prot = ""
        for line in lines:
            if line.startswith(">"):
                # Take the second token in header as the raw HLA name
                hla_names.append(line.split(" ")[1])
                if switch == 1:
                    prots.append(prot)
                    prot = ""
            else:
                prot += line.replace("\n", "")
                switch = 1
    prots.append(prot)

    hla2prot = pd.DataFrame(zip(hla_names, prots), columns=["HLA", "protin"])
    # Keep HLA like A*02:01 -> A*02:01, B*15:01 -> B*15:01
    hla2prot["HLA"] = hla2prot["HLA"].str.split(":").str[:2].str.join(":")
    hla2prot = hla2prot.drop_duplicates(subset="HLA", keep="first").reset_index(drop=True)
    return hla2prot


def aaindex_vocab(seq, vocab):
    """Encode a string into integer indices via a Vocab mapping (uppercased; 'X' -> '-')."""
    encoded = np.empty([len(seq)])
    for i in range(len(seq)):
        query = seq[i]
        if query == "X":
            query = "-"
        query = query.upper()
        encoded[i] = vocab[query]
    return encoded


class Vocab:
    """
    Minimal Vocab: maps tokens to indices and vice versa.

    Notes:
    - Index 0 is '<pad>'.
    - Unknown tokens fall back to index 0.
    """
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
            # Fall back to pad index for unknowns
            return self.token_to_idx.get(tokens, 0)
        return [self.__getitem__(token) for token in tokens]

    def to_tokens(self, indices):
        if not isinstance(indices, (list, tuple)):
            return self.idx_to_token[indices]
        return [self.idx_to_token[index] for index in indices]


# ===== Feature extraction =====
def process_chunk(chunk):
    """
    Convert a chunk of rows into feature vectors.
    Expects:
      - 'AA Composition' as a comma-separated string of floats
      - contiguous columns from 'Polarity' through 'Isoelectric Point'
    """
    chunk = chunk.reset_index(drop=True)
    features = []
    for i in range(chunk.shape[0]):
        aa_composition = chunk["AA Composition"].iloc[i].split(",")
        aa_composition = [float(value) for value in aa_composition]
        feature = aa_composition + chunk.loc[i, "Polarity":"Isoelectric Point"].tolist()
        features.append(feature)
    return features


def get_features_parallel(df, n_processes=None):
    """Parallelize feature construction across CPU processes."""
    if n_processes is None:
        n_processes = cpu_count()

    chunks = np.array_split(df, n_processes)

    with Pool(n_processes) as pool:
        results = list(tqdm(pool.imap(process_chunk, chunks), total=n_processes))

    features = [item for sublist in results for item in sublist]
    features_tensor = torch.tensor(features, dtype=torch.float32)
    return features_tensor


# ===== Parallel encoding of peptide/HLA to model inputs =====
def encode_chunk(features_chunk, peptides_chunk, hlas_chunk, hla2prot, vocab, padding_size):
    """
    Encode one chunk of examples to a fixed-length 1D vector each:
      merged = [25-dim features] + [peptide indices] + [12 zeros] + [HLA paratope indices] -> pad to padding_size
    Returns encoded list (numpy arrays) and a per-item flag (0 if HLA not found).
    """
    encoded_batch = []
    flags = []
    batch_size = len(peptides_chunk)

    for i in range(batch_size):
        peptide = peptides_chunk[i]
        hla = hlas_chunk[i]
        flag = 1

        matches = hla2prot[hla2prot["HLA"] == hla]
        if matches.shape[0] == 0:
            flag = 0
            encoded_batch.append(np.zeros(padding_size))
            flags.append(flag)
            continue

        hla_paratope = matches.iloc[0, 1] if matches.shape[0] > 0 else "-"

        peptide_encode = aaindex_vocab(peptide, vocab)
        hla_encode = aaindex_vocab(hla_paratope, vocab)
        feature_array = np.array(features_chunk[i])

        zero_array = np.zeros(12, dtype=np.float32)
        merged = np.concatenate((feature_array, peptide_encode, zero_array, hla_encode))

        current_length = merged.shape[0]
        padding_length = padding_size - current_length
        merged_padded = np.pad(merged, (0, padding_length), "constant")

        encoded_batch.append(merged_padded)
        flags.append(flag)

    return encoded_batch, flags


def encode_wrapper(args):
    return encode_chunk(*args)


def parallel_encode(features, peptides, hlas, hla2prot, vocab, padding_size, workers=4):
    """Shard inputs and encode in parallel; returns (tensor[N, L], flags)."""
    num_processes = workers
    batch_size = len(peptides)
    chunk_size = batch_size // num_processes + (batch_size % num_processes > 0)

    feature_chunks = [features[i : i + chunk_size] for i in range(0, batch_size, chunk_size)]
    peptide_chunks = [peptides[i : i + chunk_size] for i in range(0, batch_size, chunk_size)]
    hla_chunks = [hlas[i : i + chunk_size] for i in range(0, batch_size, chunk_size)]

    args_iter = (
        (fc, pc, hc, hla2prot, vocab, padding_size)
        for fc, pc, hc in zip(feature_chunks, peptide_chunks, hla_chunks)
    )

    with Pool(processes=num_processes) as pool:
        results = []
        for result in tqdm(
            pool.imap(encode_wrapper, args_iter, chunksize=max(1, len(feature_chunks) // num_processes)),
            total=len(feature_chunks),
            desc="Encoding chunks",
        ):
            results.append(result)

    encoded_arrays = [array for result in results for array in result[0]]
    flags = [flag for result in results for flag in result[1]]

    return torch.from_numpy(np.stack(encoded_arrays)), flags


# ===== Inference =====
def inference(df, thread, device, model_path, hla_prot_file, batch_size=512):
    """
    Predict immunogenicity probabilities for peptide/HLA pairs.

    Args:
      df: pandas DataFrame with columns:
          - 'peptide'
          - 'HLA'
          - 'AA Composition', 'Polarity' ... 'Isoelectric Point' (feature columns)
      thread: number of CPU threads/processes
      device: torch.device
      model_path: path to a saved model state_dict
      hla_prot_file: path to HLA paratope/protein FASTA-like file used by load_hla()
      batch_size: batch size for forward pass
    """
    padding_length = 419
    hla2prot = load_hla(hla_prot_file)

    amino = "ARNDCQEGHILKMFPSTWYV-"
    tokens = list(amino)
    vocab = Vocab(tokens, min_freq=1)

    # Limit CPU threads used by Torch ops
    torch.set_num_threads(thread)

    peptides = [p.replace(" ", "") for p in df["peptide"].astype(str).tolist()]

    def process_hla(hla):
        # Normalize when multiple alleles joined by '/', take the last segment
        hla = str(hla)
        return hla.split("/")[-1] if "/" in hla else hla

    hlas = df["HLA"].apply(process_hla).tolist()

    # Build handcrafted features
    features = get_features_parallel(df, thread)

    imm_probs = []
    print("Inference...")

    model = newBiLSTM(
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
    model = nn.DataParallel(model)

    state_dict = torch.load(model_path, map_location=try_gpu())
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        name = k.replace(".module.", ".")  # strip 'module.' if present
        new_state_dict[name] = v

    model.load_state_dict(new_state_dict)
    model.to(device=device)

    # Pre-encode entire dataset
    print("Parallel encoding...")
    encoded_all, flags_all = parallel_encode(
        features, peptides, hlas, hla2prot, vocab, padding_length, thread
    )

    # Batched forward pass
    print("Infer...")
    model.eval()
    for start in tqdm(range(0, len(peptides), batch_size)):
        end = min(start + batch_size, len(peptides))
        batch_encoded = encoded_all[start:end].to(device)
        batch_flags = flags_all[start:end]

        with torch.no_grad():
            x1 = batch_encoded[:, 25:].int()   # sequence indices region
            x2 = batch_encoded[:, 0:25].float()  # handcrafted features region
            y = model(x1, x2)
            probs = torch.softmax(y, dim=1)

            for j, flag in enumerate(batch_flags):
                imm_probs.append(float(probs[j, 1]) if flag else 0.5)

    print("Done!")
    return imm_probs


# ===== Device helpers =====
def try_gpu(i=0):
    """Return cuda:i if available; otherwise return cpu()."""
    if torch.cuda.device_count() >= i + 1:
        return torch.device(f"cuda:{i}")
    return torch.device("cpu")


def try_all_gpus():
    """Return a list of all available CUDA devices, or [cpu()] if none."""
    devices = [torch.device(f"cuda:{i}") for i in range(torch.cuda.device_count())]
    return devices if devices else [torch.device("cpu")]


# ===== HLA normalization =====
def normalize_hla(df, col="MHC Allele"):
    """
    Normalize MHC allele strings:
      1) Remove 'HLA-' prefix if present.
      2) Replace a hyphen between gene parts with a slash:
         e.g., DQA1*03:02-DQB1*04:01 -> DQA1*03:02/DQB1*04:01
      3) Re-add 'HLA-' prefix to the normalized value.
    """
    def process(allele):
        allele = str(allele).replace("HLA-", "")
        allele = re.sub(r'(?<=[0-9])-(?=[A-Z])', '/', allele)
        return f"HLA-{allele}"

    df[col] = df[col].apply(process)
    return df
