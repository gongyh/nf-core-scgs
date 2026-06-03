import numpy as np
from sklearn.cluster import KMeans
#from sklearn_extra.cluster import KMedoids
from sklearn import metrics
from sklearn.metrics import silhouette_score,davies_bouldin_score,normalized_mutual_info_score,adjusted_rand_score,fowlkes_mallows_score
from tqdm import tqdm,trange
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description="用于特征聚类:vaefeatures_dir、label_dir、fasta_dir、bins_dir")
parser.add_argument('-vd','--vaefeatures_dir',type=str,help="指定vae特征输出路径")
parser.add_argument('-ld','--label_dir',type=str,help="指定标签输出路径")
parser.add_argument('-fd','--fasta_dir',type=str,help="指定fasta文件路径")
parser.add_argument('-bd','--bins_dir',type=str,help="指定分箱结果路径")
parser.add_argument('-cvf','--cluster_num',type=str,help="分箱个数路径")
args = parser.parse_args()

with open(args.cluster_num, "r", encoding="UTF-8") as f:
    number = int(f.read().strip())
feature_path = args.vaefeatures_dir # 特征路径

# # 创建k-means模型
kmeans = KMeans(n_clusters=number,n_init=50)

 # 加载特征数据
feature_data = np.load(feature_path)

# # 对数据进行聚类
kmeans.fit(feature_data)

# # 获取每个数据点的标签
labels = kmeans.labels_
with open(args.label_dir,"w",encoding="UTF-8") as f:
    for i in labels:
        f.write(str(i)+'\n')

        
# 将聚类的结果写入到各个bin文件中
# 读取标签文件为numpy数组
labels = np.loadtxt(args.label_dir,dtype=int)
def split_data_by_label(labels, input_file, output_directory):
    bins = [[] for _ in range(max(labels) + 1)]
    # enumerate可以在每次遍历时同时获取列表中的元素和它们的索引
    # i是索引，label是元素，获取的就是对应label的索引值
    for i, label in enumerate(labels):
        bins[label].append(i)
    for i,j in enumerate(bins):
        print(f"bin{i}:{len(j)}")
    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, 'fasta'))

    for bin_num, indices in enumerate(bins):
        bin_filename = f"{output_directory}/bin{bin_num}.fasta"
        with open(bin_filename, "w") as out:
            # 将文件按列表中的索引值（列表中的元素并不是指实际的索引）写入到文件中
            for i in indices:
                out.write('>'+str(records[i].name)+'\n'+str(records[i].seq)+'\n')

input_file = args.fasta_dir
output_directory = args.bins_dir # 需要提前创建SRR952674Bins文件夹
split_data_by_label(labels, input_file, output_directory)

# ########################## 内部指标计算 ##########################
# # 计算内部指标需要传入两个文件一个是特征文件，一个是预测的标签文件
# def target(feature_Path,predLabel_Path):
#     # 加载特征数据
#     feature_data = np.load(feature_Path)

#     # 加载标签文件
#     with open(predLabel_Path ,'r',encoding='UTF-8') as f:
#         lines = f.read().splitlines()
#     silhouette_avg = silhouette_score(feature_data, lines)
#     dbi = davies_bouldin_score(feature_data,lines)

#     print("***内部指标***\n""轮廓系数：%f\nDBI指数：%f\n" %(silhouette_avg,dbi))    

# # 轮廓系数 [-1,1]值越大越好
# # DBI指数 0到正无穷，值越小越好
# f_path = "/media/ubuntu/abc/csm/myBertNew/mydata/SRR952674/usage/srr952674_unknownfeaturesxiuzheng.npy"
# l_path = "/media/ubuntu/abc/csm/myBertNew/mydata/SRR952674/usage/srr952674_unknownpredxiuzheng_num.txt"
# target(f_path,l_path)       


########################### f1,nmi等指标计算 ##########################
#################################################################

# from sklearn.metrics import silhouette_score,davies_bouldin_score,normalized_mutual_info_score,adjusted_rand_score,fowlkes_mallows_score
# #打开txt文件并读取其中的内容
# with open("/media/ubuntu/abc/csm/ablabtion_lab/sim5g/sim5g_final3_prediction.txt", 'r') as f: # 预测标签
#     content = f.read()

# with open("/media/ubuntu/abc/csm/ablabtion_lab/sim5g/sim5gover2k_real_label.txt", 'r') as fd: # 真实标签
#     contentd = fd.read()
# # 将内容按行拆分，并保存为列表
# label_pred = content.splitlines() # 预测标签
# label_true = contentd.splitlines() # 真实标签
# accuracy=metrics.accuracy_score(label_true, label_pred)
# f1=metrics.f1_score(
#           label_true, label_pred, average="macro", zero_division=0
#       )
# precision=metrics.precision_score(
#             label_true, label_pred, average="macro", zero_division=0
#         )
# recall=metrics.recall_score(
#             label_true, label_pred, average="macro", zero_division=0
#         )

# ari = adjusted_rand_score(label_true, label_pred)
# nmi = normalized_mutual_info_score(label_true, label_pred)
# fmi = fowlkes_mallows_score(label_true, label_pred)
# print("*** Pred results ***\n""acc:%f\nf1:%f\npre:%f\nrecall:%f\nari:%f\nnmi:%f\nfmi:%f\n" %(accuracy,f1,precision,recall,ari,nmi,fmi))

