#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
共现分箱核心算法
实现：覆盖度二值化 → Fisher精确检验 → t-SNE降维 → DBSCAN聚类
输入：覆盖度矩阵（TSV），super_contigs.fasta
输出：clusters.tsv（contig与簇编号），bins/（每个簇的序列）
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(description='Co-occurrence binning')
    parser.add_argument('--coverage', required=True, help='覆盖度矩阵 TSV 文件')
    parser.add_argument('--contigs', required=True, help='super_contigs.fasta 文件')
    parser.add_argument('--outdir', default='.', help='输出目录')
    parser.add_argument('--threshold', type=float, default=2048,
                        help='覆盖度二值化阈值（默认 2048 = 2^11）')
    parser.add_argument('--eps', type=float, default=2.6,
                        help='DBSCAN 邻域半径（默认 2.6）')
    parser.add_argument('--minpts', type=int, default=5,
                        help='DBSCAN 最小点数（默认 5）')
    parser.add_argument('--random_seed', type=int, default=2015,
                        help='随机种子（默认 2015）')
    return parser.parse_args()

def binarize_matrix(df, threshold):
    """将覆盖度矩阵二值化（大于阈值置1，否则0）"""
    return (df > threshold).astype(int)

def compute_pvalue_matrix(binary_mat):
    """计算所有contig对的Fisher精确检验p值矩阵（对称，对角线极小值）"""
    n = binary_mat.shape[0]
    pmat = np.ones((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            a = np.sum(binary_mat[i] & binary_mat[j])
            b = np.sum(binary_mat[i] & ~binary_mat[j])
            c = np.sum(~binary_mat[i] & binary_mat[j])
            d = np.sum(~binary_mat[i] & ~binary_mat[j])
            # Fisher exact test (two-sided)
            _, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            pmat[i, j] = p
            pmat[j, i] = p
        if (i+1) % 100 == 0:
            sys.stderr.write(f"Processed {i+1} contigs\n")
    np.fill_diagonal(pmat, 1e-14)   # 避免自距离为0
    return pmat

def tsne_d(D, random_state=2015):
    """t-SNE降维，输入距离矩阵"""
    tsne = TSNE(n_components=2, metric='precomputed', random_state=random_state)
    return tsne.fit_transform(D)

def dbscan_clustering(X, eps, min_samples):
    """DBSCAN聚类"""
    db = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean')
    return db.fit_predict(X)

def write_clusters(cluster_labels, contig_names, outdir):
    outfile = os.path.join(outdir, 'clusters.tsv')
    with open(outfile, 'w') as f:
        f.write("contig\tcluster\n")
        for name, lab in zip(contig_names, cluster_labels):
            f.write(f"{name}\t{lab}\n")
    sys.stderr.write(f"Clusters written to {outfile}\n")

def extract_bins(cluster_labels, contig_names, contig_seqs, outdir):
    bins_dir = os.path.join(outdir, 'bins')
    os.makedirs(bins_dir, exist_ok=True)
    unique_labels = set(cluster_labels)
    for lab in unique_labels:
        if lab == -1:
            continue
        outfile = os.path.join(bins_dir, f'bin_{lab}.fasta')
        with open(outfile, 'w') as f:
            for name, seq in zip(contig_names, contig_seqs):
                if cluster_labels[name] == lab:
                    f.write(f">{name}\n{seq}\n")
    sys.stderr.write(f"Bins extracted to {bins_dir}\n")

def main():
    args = parse_args()
    np.random.seed(args.random_seed)

    # 读取覆盖度矩阵
    sys.stderr.write("Loading coverage matrix...\n")
    df = pd.read_csv(args.coverage, sep='\t', index_col=0)
    contig_names = df.index.tolist()
    X = df.values
    sys.stderr.write(f"Matrix shape: {X.shape}\n")

    # 二值化
    sys.stderr.write("Binarizing...\n")
    binary = binarize_matrix(pd.DataFrame(X), args.threshold)
    sparsity = 1 - binary.sum().sum() / binary.size
    sys.stderr.write(f"Sparsity: {sparsity:.2%}\n")

    # Fisher p值矩阵
    sys.stderr.write("Computing Fisher p-values...\n")
    pmat = compute_pvalue_matrix(binary.values)

    # t-SNE
    sys.stderr.write("Running t-SNE...\n")
    X_tsne = tsne_d(pmat, random_state=args.random_seed)

    # DBSCAN
    sys.stderr.write("Running DBSCAN...\n")
    labels = dbscan_clustering(X_tsne, args.eps, args.minpts)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = np.sum(labels == -1)
    sys.stderr.write(f"Found {n_clusters} clusters, {n_noise} noise points\n")

    # 输出聚类结果
    write_clusters(labels, contig_names, args.outdir)

    # 提取bins序列
    sys.stderr.write("Extracting bin sequences...\n")
    contig_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.contigs, 'fasta')}
    seq_list = [contig_seqs[name] for name in contig_names]
    extract_bins(labels, contig_names, seq_list, args.outdir)

    sys.stderr.write("Co-occurrence binning completed.\n")

if __name__ == '__main__':
    main()
