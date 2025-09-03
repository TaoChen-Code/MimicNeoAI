# import
import collections
import pandas as pd
import re
import warnings
from collections import OrderedDict
from tqdm import tqdm

import numpy as np
from torch import nn
from torchvision import models

warnings.filterwarnings("ignore", category=FutureWarning)  # 屏蔽所有 FutureWarning


# model
class newBiLSTM(nn.Module):
    def __init__(self, vocab_size, embedding_dim, hidden_dim_x1,hidden_dim_x2, output_dim, n_layers, bidirectional, dropout, pad_idx, x2_dim):
        super(newBiLSTM, self).__init__()
        
        # 原有的嵌入和LSTM层
        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_idx)
        self.lstm = nn.LSTM(embedding_dim, hidden_dim_x1, num_layers=n_layers, bidirectional=bidirectional, dropout=dropout)
        
        # 新增的全连接层处理X2输入
        self.fc_x2 = nn.Linear(x2_dim, hidden_dim_x2)
        
        # 最终的全连接层
        self.fc = nn.Linear(hidden_dim_x1 * 2 + hidden_dim_x2, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x1, x2):
        # 在 forward 函数中调用 flatten_parameters
        self.lstm.flatten_parameters()
        # 处理X1输入
        embedded = self.dropout(self.embedding(x1.T))
        packed_output, (hidden, cell) = self.lstm(embedded)
        hidden_x1 = self.dropout(torch.cat((hidden[-2, :, :], hidden[-1, :, :]), dim=1))
        # 处理X2输入
        hidden_x2 = self.dropout(self.fc_x2(x2))
        # 拼接两个特征向量
        combined = torch.cat((hidden_x1, hidden_x2), dim=1)
        # 最终输出
        return self.fc(combined)

def load_hla(hla_prot_file):
    hla_names = []
    prots = []
    with open(hla_prot_file,'r') as f:
        lines = f.readlines()
        switch = 0
        prot = ''
        for line in lines:
            if line[0] == '>':
                hla_names.append(line.split(" ")[1])
                if switch == 1:
                    prots.append(prot)
                    prot = ''
            else:
                prot += line.replace("\n","")
                switch = 1
    prots.append(prot)
    hla2prot = pd.DataFrame(zip(hla_names,prots),columns=['HLA',"protin"])
    hla2prot['HLA'] = hla2prot['HLA'].str.split(':').str[:2].str.join(':')  # 提取前两个字段并用:连接
    hla2prot = hla2prot.drop_duplicates(subset='HLA', keep='first')  # 按HLA列去重，保留第一个出现值
    hla2prot = hla2prot.reset_index(drop=True)
    return hla2prot

def aaindex_vocab(seq,vocab):#encode peptide
    encoded = np.empty([len(seq)])#(seq_len,12)
    for i in range(len(seq)):
        query = seq[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i] = vocab[query]
    return encoded

class Vocab:
    def __init__(self, tokens=None, min_freq=0, reserved_tokens=None):
        if tokens is None:
            tokens = []
        if reserved_tokens is None:
            reserved_tokens = []
        
        # 统计词频并过滤
        counter = collections.Counter()
        for token in tokens:
            if isinstance(token, list):
                counter.update(token)
            else:
                counter[token] += 1
        
        # 按词频排序
        self.token_freqs = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        
        # 构建词汇表
        self.idx_to_token = ['<pad>'] + reserved_tokens
        self.token_to_idx = {token: idx for idx, token in enumerate(self.idx_to_token)}
        
        for token, freq in self.token_freqs:
            if freq >= min_freq and token not in self.token_to_idx:
                self.idx_to_token.append(token)
                self.token_to_idx[token] = len(self.idx_to_token) - 1
    
    def __len__(self):
        return len(self.idx_to_token)
    
    def __getitem__(self, tokens):
        if not isinstance(tokens, (list, tuple)):
            return self.token_to_idx.get(tokens, 0)  # 0对应<unk>
        return [self.__getitem__(token) for token in tokens]
    
    def to_tokens(self, indices):
        if not isinstance(indices, (list, tuple)):
            return self.idx_to_token[indices]
        return [self.idx_to_token[index] for index in indices]

## get_features_parallel
def process_chunk(chunk):
    chunk = chunk.reset_index(drop=True)
    features = []
    for i in range(chunk.shape[0]):
        aa_composition = chunk['AA Composition'].iloc[i].split(",")
        aa_composition = [float(value) for value in aa_composition]
        feature = aa_composition + chunk.loc[i, 'Polarity':'Isoelectric Point'].tolist()
        features.append(feature)
    return features

def get_features_parallel(df, n_processes=None):
    if n_processes is None:
        n_processes = cpu_count()
    
    # 将DataFrame拆分成大致相等的块
    chunks = np.array_split(df, n_processes)
    
    with Pool(n_processes) as pool:
        results = list(tqdm(pool.imap(process_chunk, chunks), total=n_processes))
    
    # 合并所有结果
    features = [item for sublist in results for item in sublist]
    features_tensor = torch.tensor(features, dtype=torch.float32)
    return features_tensor

## parallel_encode
import torch
import torch.nn.functional as F
from multiprocessing import Pool, cpu_count

def encode_chunk(features_chunk, peptides_chunk, hlas_chunk, hla2prot, vocab, padding_size):
    encoded_batch = []
    flags = []
    batch_size = len(peptides_chunk)

    for i in range(batch_size):
        peptide = peptides_chunk[i]
        hla = hlas_chunk[i]
        flag = 1
        matches = hla2prot[hla2prot['HLA'] == hla]
        if matches.shape[0] == 0:
            flag = 0
            encoded_batch.append(np.zeros(padding_size))  # 使用numpy零数组
            flags.append(flag)
            continue
        
        hla_paratope = matches.iloc[0, 1] if matches.shape[0] > 0 else "-"
        
        peptide_encode = aaindex_vocab(peptide, vocab)  # 直接返回numpy数组
        hla_encode = aaindex_vocab(hla_paratope, vocab)
        feature_array = np.array(features_chunk[i])
        
        # 使用numpy进行拼接
        zero_array = np.zeros(12, dtype=np.float32)
        merged = np.concatenate((feature_array, peptide_encode, zero_array, hla_encode))

        current_length = merged.shape[0]
        padding_length = padding_size - current_length
        merged_padded = np.pad(merged, (0, padding_length), 'constant')

        encoded_batch.append(merged_padded)
        flags.append(flag)

    return encoded_batch, flags

def encode_wrapper(args):
    return encode_chunk(*args)

def parallel_encode(features, peptides, hlas, hla2prot, vocab, padding_size, workers=4):
    num_processes = workers
    batch_size = len(peptides)
    chunk_size = batch_size // num_processes + (batch_size % num_processes > 0)

    # 数据分块
    feature_chunks = [features[i:i + chunk_size] for i in range(0, batch_size, chunk_size)]
    peptide_chunks = [peptides[i:i + chunk_size] for i in range(0, batch_size, chunk_size)]
    hla_chunks = [hlas[i:i + chunk_size] for i in range(0, batch_size, chunk_size)]

    # 生成参数迭代器
    args_iter = ((fc, pc, hc, hla2prot, vocab, padding_size) 
                for fc, pc, hc in zip(feature_chunks, peptide_chunks, hla_chunks))

    with Pool(processes=num_processes) as pool:
        # 使用imap+进度条
        results = []
        for result in tqdm(pool.imap(
            encode_wrapper,
            args_iter,
            chunksize=max(1, len(feature_chunks)//num_processes)
        ), total=len(feature_chunks), desc="Encoding chunks"):
            results.append(result)

    # 合并结果
    encoded_arrays = [array for result in results for array in result[0]]
    flags = [flag for result in results for flag in result[1]]
    
    return torch.from_numpy(np.stack(encoded_arrays)), flags


def inference(df, thread, device, model_path, batch_size=512):
    padding_length = 419
    hla2prot = load_hla()
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    tokens = list(amino)
    vocab = Vocab(tokens, min_freq=1)

    # 设置 CPU 使用的线程数量
    torch.set_num_threads(thread)
    peptides = df['peptide'].tolist()
    peptides = [p.replace(" ", "") for p in peptides]

    def process_hla(hla):
        if '/' in hla:
            return hla.split('/')[-1]
        else:
            return hla

    hlas = df['HLA'].apply(process_hla).tolist()

    # 获取特征
    features = get_features_parallel(df, thread)

    imm_probs = []
    print("inference....")

    model = newBiLSTM(vocab_size=len(vocab), embedding_dim=12, hidden_dim_x1=256,
                      hidden_dim_x2=16, output_dim=2, n_layers=2, bidirectional=True,
                      dropout=0.5, pad_idx=vocab['<pad>'], x2_dim=25)
    model = nn.DataParallel(model)

    state_dict = torch.load(model_path, map_location=try_gpu())
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        name = k.replace('.module.', '.')  # 移除 `module.`
        new_state_dict[name] = v

    model.load_state_dict(new_state_dict)
    model.to(device=try_gpu())

    # 预处理阶段：全量数据并行编码
    print("Parallel encoding...")
    encoded_all, flags_all = parallel_encode(
        features, peptides, hlas, hla2prot, vocab, padding_length, thread
    )

    # 按照 batch_size 分批处理
    print("Infer...")
    model.eval()
    for start in tqdm(range(0, len(peptides), batch_size)):
        end = min(start + batch_size, len(peptides))
        batch_encoded = encoded_all[start:end].to(device)
        batch_flags = flags_all[start:end]

        with torch.no_grad():
            x1 = batch_encoded[:, 25:].int()
            x2 = batch_encoded[:, 0:25].float()
            y = model(x1, x2)
            probs = torch.softmax(y, dim=1)

            for j, flag in enumerate(batch_flags):
                imm_probs.append(float(probs[j, 1]) if flag else 0.5)

    print("done!")
    return imm_probs


def try_gpu(i=0):  #@save
    """如果存在，则返回gpu(i)，否则返回cpu()"""
    if torch.cuda.device_count() >= i + 1:
        return torch.device(f'cuda:{i}')
    return torch.device('cpu')

def try_all_gpus():  #@save
    """返回所有可用的GPU，如果没有GPU，则返回[cpu(),]"""
    devices = [torch.device(f'cuda:{i}')
             for i in range(torch.cuda.device_count())]
    return devices if devices else [torch.device('cpu')]



def normalize_hla(df, col="MHC Allele"):
    """
    处理 MHC Allele 列：
    1. 去掉前缀 HLA-
    2. 用正则将基因片段之间的 '-' 改成 '/'
       例如: DQA1*03:02-DQB1*04:01 -> DQA1*03:02/DQB1*04:01
    3. 给所有结果重新添加前缀 HLA-
    """
    def process(allele):
        allele = str(allele).replace("HLA-", "")  # 去掉已有前缀
        # 用正则替换：在等位基因片段之间的 '-' 改为 '/'
        allele = re.sub(r'(?<=[0-9])-(?=[A-Z])', '/', allele)
        return f"HLA-{allele}"

    df[col] = df[col].apply(process)
    return df


