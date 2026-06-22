#!/usr/bin/env python3

import os
import pathlib
import logging
from sklearn.cluster import KMeans
import numpy as np
# from Bio import SeqIO
import logging
import os
import argparse
import sys
import math
import re
import csv
import time
from sklearn.preprocessing import normalize

"""
文件运行路径需要在marker_gene下
且需要切换到人院服务器下的test环境下
yanziming@server1:~/vicent/contrast_experiment/marker_gene$ 
"""
# create logger
logger = logging.getLogger('ourmethod')


def silhouette(X, W, label):
    X_colsum = np.sum(X ** 2, axis=1)
    X_colsum = X_colsum.reshape(len(X_colsum), 1)
    W_colsum = np.sum(W ** 2, axis=1)
    W_colsum = W_colsum.reshape(len(W_colsum), 1)

    Dsquare = np.tile(
        X_colsum, (1, W.shape[0])) + np.tile(W_colsum.T, (X.shape[0], 1)) - 2 * X.dot(W.T)
    # avoid error caused by accuracy
    Dsquare[Dsquare < 0] = 0
    D = np.sqrt(Dsquare)
    aArr = D[np.arange(D.shape[0]), label]
    D[np.arange(D.shape[0]), label] = np.inf
    bArr = np.min(D, axis=1)
    tmp = (bArr - aArr) / np.maximum(aArr, bArr)
    return np.mean(tmp)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# Modified from SolidBin


def scan_for_marker_genes(contig_file, x_contigs, nthreads, bestK=0):
    # 获取路径
    software_path = pathlib.Path(__file__).parent.parent.absolute()
    fragScanURL = os.path.join(software_path, 'auxiliary',
                               'FragGeneScan1.31', 'run_FragGeneScan.pl')
    hmmExeURL = os.path.join(software_path, 'auxiliary', 'hmmer-3.3',
                             'src', 'hmmsearch')
    markerExeURL = os.path.join(
        software_path, 'auxiliary', 'test_getmarker.pl')
    markerURL = os.path.join(software_path, 'auxiliary', 'marker.hmm')
    logger.debug(markerURL)

    seedURL = contig_file + ".seed"
    fragResultURL = contig_file+".frag.faa"
    hmmResultURL = contig_file+".hmmout"
    
    candK = 3
    maxK = min(20, len(x_contigs))
    stepK = 2
    
    # 文件不在就使用脚本生成
    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL+" -genome="+contig_file+" -out="+contig_file + \
            ".frag -complete=0 -train=complete -thread="+str(nthreads)+" 1>" + \
            contig_file+".frag.out 2>"+contig_file+".frag.err"
        logger.debug("exec cmd: "+fragCmd)
        os.system(fragCmd)
    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL+" --domtblout "+hmmResultURL+" --cut_tc --cpu "+str(nthreads)+" " + \
                markerURL+" "+fragResultURL+" 1>"+hmmResultURL+".out 2>"+hmmResultURL+".err"
            logger.debug("exec cmd: "+hmmCmd)
            os.system(hmmCmd)
        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + \
                    " " + contig_file + " 1000 " + seedURL
                logger.debug("exec cmd: " + markerCmd)
                os.system(markerCmd)
            if os.path.exists(seedURL):
                candK = file_len(seedURL)+2
                maxK = min(3 * candK, len(x_contigs))
                stepK = 2

            else:
                logger.info("seed not exist, k start from 3 ")
                candK = 3
                maxK = min(20, len(x_contigs))
                stepK = 2

        else:
            logger.debug("HMMER search failed! Path: " +
                         hmmResultURL + " does not exist.")
    else:
        logger.debug("FragGeneScan failed! Path: " +
                     fragResultURL + " does not exist.")
    # X_mat=x_contigs.numpy()
    X_mat = x_contigs
    if bestK == 0:
        bestK = candK
        if candK == maxK:
            bestK = candK
            bestSilVal = 0
        else:
            bestSilVal = 0
            t = time.time()
            for k in range(candK, maxK, stepK):
                kmeans = KMeans(n_clusters=k, init='k-means++',
                                random_state=9, n_jobs=-1)
                kmeans.fit(X_mat)
                silVal = silhouette(
                    X_mat, kmeans.cluster_centers_, kmeans.labels_)
                logger.info("k:" + str(k) + "\tsilhouette:" +
                            str(silVal) + "\telapsed time:" + str(time.time() - t))
                t = time.time()

                if silVal > bestSilVal:
                    bestSilVal = silVal
                    bestK = k
                else:
                    break

        # candKold=candK
        candK = bestK + candK
        if maxK > candK:
            bestSilVal_2nd = 0
            for k in range(candK, maxK, stepK):
                kmeans = KMeans(n_clusters=k, init='k-means++',
                                random_state=9, n_jobs=-1)
                kmeans.fit(X_mat)
                silVal_2nd = silhouette(
                    X_mat, kmeans.cluster_centers_, kmeans.labels_)
                logger.info("k:" + str(k) + "\tsilhouette:" +
                            str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
                t = time.time()
                if silVal_2nd > bestSilVal_2nd:
                    bestSilVal_2nd = silVal_2nd
                    bestK_2nd = k
                else:
                    break
            if bestSilVal_2nd > bestSilVal:
                bestSilVal = bestSilVal_2nd
                bestK = bestK_2nd

        logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))

    return bestK

# Get contigs containing marker genes


def get_contigs_with_marker_genes(mg_length_threshold,
                                  contig_lengths,
                                  min_length):

    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(f"/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/shotgun/SRR1262938/spades_out/contigs_over_1000.fasta.hmmout", "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                # Get contig name
                name_strings = contig.split("_")
                name_strings = name_strings[:len(name_strings)-3]
                contig_name = "_".join(name_strings)

                contig_length = contig_lengths[contig_name]

                if contig_length >= min_length and mapped_marker_length > marker_gene_length*mg_length_threshold:

                    marker_repeated_in_contig = False

                    # Get marker genes in each contig
                    if contig_name not in contig_markers:
                        contig_markers[contig_name] = [marker_gene]
                    else:
                        if marker_gene not in contig_markers[contig_name]:
                            contig_markers[contig_name].append(marker_gene)

                    # Get contigs containing each marker gene
                    if marker_gene not in marker_contigs:
                        marker_contigs[marker_gene] = [contig_name]
                    else:
                        if contig_name not in marker_contigs[marker_gene]:
                            marker_contigs[marker_gene].append(contig_name)
                        else:
                            marker_repeated_in_contig = True

                    # Get contig counts for each marker
                    if marker_gene not in marker_contig_counts:
                        marker_contig_counts[marker_gene] = 1
                    else:
                        if not marker_repeated_in_contig:
                            marker_contig_counts[marker_gene] += 1

    return marker_contigs, marker_contig_counts, contig_markers


def count_contigs_with_marker_genes(marker_contig_counts):

    marker_frequencies = {}

    for marker in marker_contig_counts:

        if marker_contig_counts[marker] not in marker_frequencies:
            marker_frequencies[marker_contig_counts[marker]] = 1
        else:
            marker_frequencies[marker_contig_counts[marker]] += 1

    return marker_frequencies


def pairwise_combinations(input_list):
    return [(input_list[i], input_list[j]) for i in range(len(input_list)) for j in range(i + 1, len(input_list))]


def write_2d_list_to_txt(data_2d_list, file_path):
    with open(file_path, 'w') as file:
        for row in data_2d_list:
            line = ' '.join(str(element) for element in row)
            file.write(line + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-kf', '--kmer_file', help='kmer路径')
    parser.add_argument('-cf', '--configs_file',help='configs_file路径')
    parser.add_argument('-cvf', '--cluster_value',help='聚类个数输出路径')
    args = parser.parse_args()
    # 节点特征文件
    # x_contigs = np.loadtxt("/home/yanziming/vicent/SGCMCV2/2021-TMM-SGCMC/data/所有contigs的特征.txt",dtype=float)
    # 读取csv文件
    x_contigs = np.loadtxt(args.kmer_file,
                           delimiter=',', usecols=range(1, 137))

    # x_contigs = normalize(x_contigs, norm='l2', axis=1) # 会影响最终的结果，全部序列归一化后结果为12，归一化之前是14。未知序列归一化后结果为4，归一化之前也为4。
    # 保存数据用来欧拉变换
    # np.savetxt(
    #     '/home/yanziming/csm/data_set/sharon/X1.txt', x_contigs)
    # fasta文件
    # contig_file = '/home/yanziming/vicent/data_set/human_gut_metagenome/shotgun/SRR6131123/spades_out/contigs_over_1000.fasta'
    
    # 最终确定的聚类个数
    cluster_k = scan_for_marker_genes(
        args.configs_file, x_contigs, 40, 0)  # 1为脚本运行线程数
    with open(args.cluster_value, "w", encoding="UTF-8") as f:
        value = cluster_k * 1.5
        int_value = int(math.ceil(value))  # 向上取整再取整数部分
        unit = int_value % 10
        if unit <= 5:
            rounded_value = int_value - unit + 5
        else:
            rounded_value = int_value - unit + 10
        f.write(str(rounded_value))
    
    

# # # 提取数量
#     contig_length = {}
#     for record in SeqIO.parse(contig_file, "fasta"):
#         contig_length[record.id] = len(record.seq)
#     # import pickle
#     # pickle.dump(contig_length, open(
#     #     f"{output}/profiles/contig_lengths.pkl", "wb+")) #改文件路径
#     # scan_for_marker_genes(contigs, output, threads) #contigs，输出，线程
#     marker_contigs, marker_contig_counts, contig_markers = get_contigs_with_marker_genes(
#         "./yeast", 0.5, contig_length, 1000)
#     print(type(marker_contigs))
#     print(len(marker_contigs))
#     print(type(marker_contig_counts))
#     print(type(contig_markers))
#     # file1 = open("./yeast/marker_contigs", 'w')
#     # file1.write(str(marker_contigs))
#     # file2 = open("./yeast/marker_contig_counts", 'w')
#     # file2.write(str(marker_contig_counts))
#     # file1 = open("./yeast/contig_markers", 'w')
#     # file1.write(str(contig_markers))

#     # 处理marker_contigs字典为约束型M//处理M约束完成
#     # 处理每个键值对
#     node_set_list = []
#     id = []
#     id_list = []
#     for key, values in marker_contigs.items():
#         # print(type(values))
#         # break
#         id = []
#         for i in range(len(values)):
#             tmp = values[i].split('_')
#             values[i] = '_'.join(tmp[:2])
#             id.append(tmp[1])  # 存储数字
#         id_list.append(id)  # 存储节点数字列表
#     n = 0
# # # 示例用法
# # my_list = [1, 2, 3, 4]
# # combinations = pairwise_combinations(my_list)

# # # 将所有组合放入一个二维列表
# # result_2d_list = [list(pair) for pair in combinations]

# # # 打印结果
# # print(result_2d_list)
#     # 遍历id_list列表两两组合找到互斥连接
#     must_not_connect = []
#     for i in range(len(id_list)):
#         if(len(id_list[i]) == 1):
#             continue
#         combinations = pairwise_combinations(id_list[i])
#         combinations = [list(pair) for pair in combinations]
#         for i in range(len(combinations)):
#             must_not_connect.append(combinations[i])
#     must_not_connect2 = []
#     for i in range(len(must_not_connect)):
#         must_not_connect2.append(must_not_connect[i][::-1])
#     must_not_connect = must_not_connect+must_not_connect2
#     print(must_not_connect[2])
#     print(f"互斥连接的个数为{len(must_not_connect)}")
#     # data_2d_list = [[1, 2], [3, 4], [5, 6]]
#     file_path = '/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/hic_map/must_not_contact.txt'
#     write_2d_list_to_txt(must_not_connect, file_path)
#     # /home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/hic_map/hic_map_without_weight.txt
#     # 打开文件并逐行读取内容
#     with open('/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/hic_map/hic_map_without_weight.txt', 'r') as file:
#         lines = file.readlines()

#     # 打印逐行读取到的内容
#     maybe_contact = []
#     for line in lines:
#         # print(line.strip())  # 使用strip()方法移除行尾的换行符
#         maybe_contact.append(line.strip().split())
#     # print(maybe_contact[0])
#     # print(type(maybe_contact[0]))

#     result = [
#         sublist for sublist in maybe_contact if sublist not in must_not_connect]
#     print(len(result))
#     print(len(must_not_connect))
#     print(len(maybe_contact))
#     file_path = '/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/hic_map/must_contact.txt'
#     write_2d_list_to_txt(result, file_path)
#     # # 举例说明
#     # list1 = [['1', '2', '3'], ['4', '5', '6'],
#     #          ['7', '8', '9'], ['11', '12', '13']]
#     # list2 = [['1', '2', '3'], ['4', '5', '6'], ['10', '11', '12']]

#     # # 使用列表推导式查找list1中存在而list2中不存在的子列表
#     # result = [sublist for sublist in list1 if sublist not in list2]

#     # print(result)  # 输出: [[7, 8, 9], [11, 12, 13]]

#     for key, values in marker_contigs.items():
#         tmp = ', '.join(values)
#         node_set_list.append('Set '+str(n)+': '+tmp)
#         n += 1
#     # print(id_list)
#     # markers_list = ['marker1', 'marker2', 'marker3']

#     # with open('./yeast/contigs.fasta.markers', 'w') as file:
#     #     for node_set in node_set_list:
#     #         file.write(node_set + '\n')

#     # print(marker_contigs)
