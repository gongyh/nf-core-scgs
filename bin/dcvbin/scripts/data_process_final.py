# 去掉final.contigs.fa中2kbp以下的序列###########
from Bio import SeqIO
from tqdm import tqdm
import itertools
import argparse
import random
import csv

parser = argparse.ArgumentParser(description="用于数据集预处理，截取和构造进一步预训练样本，保留序列")
parser.add_argument('-fasta','--initial_fasta',type=str,help="原始fasta路径")
parser.add_argument('-cfasta','--cut_fasta',type=str,help="截取2k后的fasta路径")
parser.add_argument('-aseq','--all_seq',type=str,help="只保留序列用于特征提取")
parser.add_argument('-ftseq','--futher_train_seq',type=str,help="提取进一步预训练序列")
parser.add_argument('-tseq','--train_seq',type=str,help="进一步预训练训练数据集")
parser.add_argument('-dseq','--dev_seq',type=str,help="进一步预训练验证数据集")
args = parser.parse_args()

with open(args.fasta, 'r') as infile:
    with open(args.cfasta, 'w') as outfile:
    # 使用SeqIO.parse读取FASTA文件中的序列
        for record in SeqIO.parse(infile, 'fasta'):
            # 检查序列长度是否大于或等于阈值
            if len(record.seq) >= 2000:
                # 如果满足条件，将序列写入输出文件
                SeqIO.write(record, outfile, 'fasta')

with open(args.aseq,"w") as ouput_file:
    for record in tqdm(SeqIO.parse(args.cfasta,"fasta")):
        seq = str(record.seq)
        ouput_file.write(seq+'\n')

# 读取FASTA文件序列信息并存入列表
records = [str(record.seq) for record in SeqIO.parse(args.cfasta, "fasta")]
with open(args.ftseq,"w",newline = '') as csvfile:
    writer = csv.writer(csvfile)
    for sequence in tqdm(records,desc = "正在生成数据集："):
        segments = [sequence[i:i+10000] for i in range(0, len(sequence), 10000) if len(sequence[i:i+10000]) == 10000]
        # 生成两两组合
        combinations = list(itertools.combinations(segments, 2))
        writer.writerows(combinations)

# 拆分数据集
# 从CSV文件中读取所有行
all_rows = []
with open(args.ftseq, "r", newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        all_rows.append(row)

# 随机打乱行的顺序
random.shuffle(all_rows)

# 根据8:2的比例计算分割索引
split_index = int(0.8 * len(all_rows))

# 分割为训练集和测试集
train_set = all_rows[:split_index]
val_set = all_rows[split_index:]

# 将训练集写入CSV文件
with open(args.tseq, "w", newline='') as train_csv:
    writer = csv.writer(train_csv)
    writer.writerows(train_set)

# 将测试集写入CSV文件
with open(args.dseq, "w", newline='') as val_csv:
    writer = csv.writer(val_csv)
    writer.writerows(val_set)


# 检测可以生成多少条数
def calculate_combinations(n):
    return n * (n - 1) // 2

total_combinations = 0

# 使用tqdm显示进度条
for record in tqdm(SeqIO.parse(args.cfasta, "fasta")):
    seq = str(record.seq)
    if len(seq) >= 20000:
        i = len(seq) // 10000
        s = calculate_combinations(i)
        #print(f"Sequence length: {len(seq)}, Segments: {i}, Combinations: {s}")
        total_combinations += s

print(f"Total combinations: {total_combinations}")