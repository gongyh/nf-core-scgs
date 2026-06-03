import argparse
from Bio import SeqIO
from tqdm import tqdm
import itertools
import random
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument('-fasta', type=str, required=True, help="输入的fasta文件路径")
parser.add_argument('-outdir', type=str, required=True, help="输出目录")
args = parser.parse_args()

# 读取序列
records = [str(record.seq) for record in SeqIO.parse(args.fasta, "fasta")]

# 写入所有组合对
fur_output = os.path.join(args.outdir, "furtrain_output.csv")
with open(fur_output, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    for sequence in tqdm(records, desc="正在生成预训练数据对："):
        segments = [sequence[i:i+20000] for i in range(0, len(sequence), 20000) if len(sequence[i:i+20000]) == 20000]
        combinations = list(itertools.combinations(segments, 2))
        writer.writerows(combinations)

# 拆分数据集
with open(fur_output, "r", newline='') as csvfile:
    all_rows = list(csv.reader(csvfile))
random.shuffle(all_rows)
split_index = int(0.9 * len(all_rows))
train_set = all_rows[:split_index]
val_set = all_rows[split_index:]

with open(os.path.join(args.outdir, "train_set.csv"), "w", newline='') as train_csv:
    csv.writer(train_csv).writerows(train_set)
with open(os.path.join(args.outdir, "val_set.csv"), "w", newline='') as val_csv:
    csv.writer(val_csv).writerows(val_set)
