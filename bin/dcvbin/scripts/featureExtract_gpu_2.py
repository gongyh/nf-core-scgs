#SequentialSampler# import torch
# import time
# import numpy as np
# import argparse
# from tqdm import tqdm
# from Bio import SeqIO
# from transformers import AutoTokenizer,AutoModel
import torch
import os
import time
import pandas as pd
import numpy as np
import random
from sklearn.metrics import f1_score
from Bio import SeqIO
from transformers import BertModel, BertTokenizer,BertConfig,BertForPreTraining 
from transformers import BertForSequenceClassification as bert
from torch.utils.data import DataLoader, RandomSampler, SequentialSampler
from transformers import AdamW, get_linear_schedule_with_warmup
from sklearn import metrics
from tqdm import tqdm
from copy import deepcopy
from sklearn.model_selection import train_test_split
from torch.utils.data import TensorDataset # 处理数据集
from transformers import AutoTokenizer,AutoModel
import time
import argparse
# 设置参数
parser = argparse.ArgumentParser(description="用于特征提取:model_dir、seq.txt、feature_output.npy")
parser.add_argument('-md','--model_dir',type=str,help="指定模型路径")
parser.add_argument('-fd','--fasta_file',type=str,help="指定输入序列路径")
parser.add_argument('-sd','--seq_file',type=str,help="提取序列路径")
parser.add_argument('-dd','--fpf_file',type=str,help="指定特征输出路径")
args = parser.parse_args()

# 设置设备为GPU（如果可用）
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


#tokenizer = AutoTokenizer.from_pretrained(args.model_dir, trust_remote_code=True)
#model = AutoModel.from_pretrained(args.model_dir, trust_remote_code=True).to(device)

from pathlib import Path

# 统一解析模型目录为绝对路径（支持相对路径）
model_dir = Path(args.model_dir).resolve()

#tokenizer = AutoTokenizer.from_pretrained(str(model_dir), trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained(str(model_dir), trust_remote_code=True, local_files_only=True)
model = AutoModel.from_pretrained(str(model_dir), trust_remote_code=True).to(device)


with open(args.seq_file,"w") as ouput_file:
    for record in tqdm(SeqIO.parse(args.fasta_file,"fasta")):
        seq = str(record.seq)
        if len(seq)>=2000:
            ouput_file.write(seq+'\n')

time.sleep(5)

# Read sequences from file#
with open(args.seq_file, "r", encoding="UTF-8") as f:
    seq_lines = f.read().splitlines()

feature_list = []
start_time = time.time()
for seq in tqdm(seq_lines,desc = "正在提取特征......"):
    inputs = tokenizer(seq, 
                       return_tensors = 'pt',
                       padding="longest",
                       max_length=5000,
            truncation=True,)["input_ids"]
    inputs = inputs.to(device)  # 将输入数据移动到GPU
    with torch.no_grad():  # 禁用梯度计算以节省内存
        hidden_states = model(inputs)[0] # [1, sequence_length, 768]
        # embedding with mean pooling
        embedding_mean = torch.mean(hidden_states[0], dim=0)
        feature_list.append(embedding_mean.detach().cpu().numpy())

np.save(args.fpf_file, np.stack(feature_list))
end_time = time.time()
print("提取特征耗费时间：%d" %(end_time-start_time))
print(len(feature_list))





