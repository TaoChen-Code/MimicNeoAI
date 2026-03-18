import collections
from collections import OrderedDict

import torch
from torch import nn


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


class NewBiLSTM(nn.Module):
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
        super().__init__()

        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_idx)
        self.lstm = nn.LSTM(
            embedding_dim,
            hidden_dim_x1,
            num_layers=n_layers,
            bidirectional=bidirectional,
            dropout=dropout,
        )
        self.fc_x2 = nn.Linear(x2_dim, hidden_dim_x2)
        self.fc = nn.Linear(hidden_dim_x1 * 2 + hidden_dim_x2, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x1, x2):
        self.lstm.flatten_parameters()
        embedded = self.dropout(self.embedding(x1.T))
        _, (hidden, _) = self.lstm(embedded)
        hidden_x1 = self.dropout(torch.cat((hidden[-2, :, :], hidden[-1, :, :]), dim=1))
        hidden_x2 = self.dropout(self.fc_x2(x2))
        combined = torch.cat((hidden_x1, hidden_x2), dim=1)
        return self.fc(combined)



def build_model(vocab_size: int, pad_idx: int) -> nn.Module:
    return NewBiLSTM(
        vocab_size=vocab_size,
        embedding_dim=12,
        hidden_dim_x1=256,
        hidden_dim_x2=16,
        output_dim=2,
        n_layers=2,
        bidirectional=True,
        dropout=0.5,
        pad_idx=pad_idx,
        x2_dim=25,
    )



def load_model_weights(model: nn.Module, model_path: str, device: torch.device) -> None:
    state_dict = torch.load(model_path, map_location=device)

    candidate_dicts = [state_dict]

    stripped = OrderedDict()
    for k, v in state_dict.items():
        stripped[k.replace("module.", "", 1)] = v
    candidate_dicts.append(stripped)

    prefixed = OrderedDict()
    for k, v in state_dict.items():
        if not k.startswith("module."):
            prefixed[f"module.{k}"] = v
        else:
            prefixed[k] = v
    candidate_dicts.append(prefixed)

    last_error = None
    for cand in candidate_dicts:
        try:
            model.load_state_dict(cand, strict=True)
            return
        except RuntimeError as exc:
            last_error = exc

    raise RuntimeError(
        f"Failed to load model weights from {model_path}. Last error: {last_error}"
    )
