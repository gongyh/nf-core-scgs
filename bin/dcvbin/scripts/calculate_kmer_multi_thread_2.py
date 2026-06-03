import threading
import pandas as pd
import numpy as np
import getopt
import sys
from tqdm import tqdm
import csv
import argparse



"""
此脚本的使用方法：
python calculate_kmer_multi_thread.py <input_fasta_path> <kmer_output_path> <seqId_output_path> <k_value>
参数说明：原始fasta文件，kmer频率存储路径，提取的序列iD存储路径，设置的kmer频率
例子：python calculate_kmer_multi_thread.py "unknownseq.fasta" "result\kmer.csv" "result\SeqId.csv" 4
python calculate_kmer_multi_thread.py SRR17858159over2k.fasta 4mer.csv SeqId.csv 4
cmd中直接粘贴路径即可，不用对'\4'或'\t'之类的进行转义。
"""

def caculate_kmer(k):
    #opts, args = getopt.getopt(sys.argv[1:], "hi:k:")
    all_title = []
    all_contig = []
    # %(vaa,avv)
    # input_file = "/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/shotgun/SRR1262938/megahit_out/final.contigs.fa"
    # intput_file = args.input
    #        input_file="D:\\宏基因组101\\序列\\0+0.txt"
    # 选择kmer的长度
    # k = 4
    # 输入文件
    with open(args.input) as f:
        list = f.readlines()

    N = 0  # number of contigs
    # split contigs title and sequence into different lists

    a = ""
    #        title=[]
    contig = []
    #    title.append(list[0])
    # N+=1
    for i in tqdm(range(int(len(list)))):
        list[i] = list[i].strip('\n')
        if (">" in list[i]):
            all_title.append(list[i])
            if a != "":
                contig.append(a)
            a = ""
            N += 1
        else:
            a += list[i]
    contig.append(a)
    #        all_title.append(title)
    seq = []
    for i in tqdm(range(0, len(contig))):
        x = str(contig[i]).replace("A", "0").replace("a", "0").replace("T", "1").replace("t", "1").replace("G", "2").replace("g", "2").replace("C", "3").replace("c", "3").replace(",", "").replace("F", "0").replace(
            "R", "2").replace("Y", "3").replace("M", "0").replace("K", "1").replace("S", "3").replace("W", "0").replace("H", "1").replace("B", "1").replace("D", "2").replace("V", "0").replace("N", "0")
        seq.append(x)

    vec = np.zeros(pow(4, k))  # k^4，初始化vec向量
    # mer=get_list(line,4) # k=4，获取k-mer片段存入mer
    weight = []
    for i in range(k):
        weight.append(pow(4, i))  # 4^(k-1),4^(k-2)……4^0，设定权重向量
    weight = weight[::-1]  # [64, 16, 4, 1]
    for j in tqdm(range(len(seq))):
        # mer[j]
        # mm=list(mer[j]) # 将一串字符串打散
        m = seq[j]
        for i in range(0, len(seq[j])-k+1):  # 将字符串中的数字分别转换成整数并形成矩阵
            m1 = m[i:i+k]
            a = 0
            for ii in range(len(m1)):
                if m1[ii] not in ["0", "1", '2', '3', '4', '5', '6']:
                    print(m1)
                a += weight[ii]*int(m1[ii])
            vec[a] = vec[a]+1  # 得到256维的k-mer频率向量vec
        # 反向互补
    #print (vec)

        for ii in range(len(vec)):
            if vec[ii] == -1:
                continue
            else:
                '''
                a=int(ii/64)
                b=int((ii-a*64)/16)
                c=int((ii-a*64-b*16)/4)
                d=ii-a*64-b*16-c*4 
                rank=np.array([a,b,c,d]) # 转化成4位二进制数
                '''
                ind = ii
                rank = np.zeros(k)
                for jj in range(k):
                    rank[jj] = int(ind/pow(4, k-jj-1))
                    ind -= rank[jj]*pow(4, k-jj-1)

                rank1 = np.zeros(k)
                aa = np.argwhere(rank == 0)  # 换成互补序列
                bb = np.argwhere(rank == 1)
                cc = np.argwhere(rank == 2)
                dd = np.argwhere(rank == 3)
                for jj in range(len(aa)):
                    rank1[aa[jj]] = 1
                for jj in range(len(bb)):
                    rank1[bb[jj]] = 0
                for jj in range(len(cc)):
                    rank1[cc[jj]] = 3
                for jj in range(len(dd)):
                    rank1[dd[jj]] = 2
                rank1 = rank1[::-1]  # 倒序输出
                i1 = 0
                for jj in range(len(rank1)):
                    i1 += rank1[jj]*weight[jj]
                if int(i1) != ii:
                    ee = vec[ii]+vec[int(i1)]
                    vec[ii] = ee
                    vec[int(i1)] = -1
        vecf = vec[vec != -1]  # 得到最后的k-mer频率矩阵vecf，共136维
        vec = np.zeros(pow(4, k))
    # vecm=np.row_stack((vecm,vecf))
    #            X.append(vecf)
        all_contig.append(vecf)

    # X1=np.array(X,dtype=float)
    #kmer = [[float(x) for x in y.split(',') if len(x) >= 1] for y in X1[1:] if len(y) >= 1]
    test = pd.DataFrame(data=all_contig)
    k = str(k)
    # out_file = '/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/process_data/4ker.csv'
    test.to_csv(args.outputkmer, encoding='gbk', header=0)

    test2 = pd.DataFrame(data=all_title)
    # out_file2 = '/home/yanziming/vicent/data_set/synthetic_metagenomic_yeast/process_data/4ker_title.csv'
    test2.to_csv(args.outputkmertitle, encoding='gbk', header=0)

    print("Successful! Please go csv for a look.")


if __name__ == '__main__':
    # 初始化参数构造器
    parser = argparse.ArgumentParser()
    # 添加参数
    parser.add_argument('input', help='input.fasta路径')
    parser.add_argument('outputkmer', help='kmer.csv路径')
    parser.add_argument('outputkmertitle',help='seqId.csv路径')
    parser.add_argument('knum', help='kmer长度')

    # 解析参数
    args = parser.parse_args()
    t = threading.Thread(target=caculate_kmer(int(args.knum)))
    t.start()  # 启动线程，让线程开始执行
    t.join()  # 等待线程执行完毕
