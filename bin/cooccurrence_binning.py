#!/usr/bin/env python3
"""
Co-occurrence binning based on Fisher's exact test and DBSCAN.
Input: coverage matrix (TSV: contigs × samples, values = coverage)
Output: clusters.tsv (contig -> bin)
Usage: python cooccurrence_binning.py coverage.tsv clusters.tsv [--eps EPS] [--min_samples MIN_SAMPLES]
"""

import sys
import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(description='Co-occurrence binning')
    parser.add_argument('coverage_file', help='coverage matrix TSV')
    parser.add_argument('filtered_fasta', help='filtered contigs FASTA file')
    parser.add_argument('output_file', help='output clusters TSV')
    parser.add_argument('--eps', type=float, default=0.05,
                        help='DBSCAN eps (default: 0.05)')
    parser.add_argument('--min_samples', type=int, default=2,
                        help='DBSCAN min_samples (default: 2)')
    parser.add_argument('--tsne', action='store_true',
                        help='Apply t-SNE before DBSCAN')
    parser.add_argument('--tsne_dim', type=int, default=3,
                        help='t-SNE output dimensions (default: 3)')
    parser.add_argument('--tsne_perplexity', type=float, default=30.0,
                        help='t-SNE perplexity (default: 30)')
    args = parser.parse_args()
    contig_ids = []
    with open(args.filtered_fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line[1:].split()[0]
                contig_ids.append(contig_id)
    print(f"Filtered contigs from fasta: {len(contig_ids)}", file=sys.stderr)
    df = pd.read_csv(args.coverage_file, sep='\t', index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    valid_ids = [id for id in contig_ids if id in df.index]
    print(f"Valid contigs in matrix: {len(valid_ids)}", file=sys.stderr)
    if len(valid_ids) == 0:
        print("No valid contigs found after filtering.", file=sys.stderr)
        with open(args.output_file, 'w') as f:
            f.write("contig\tbin\n")
        return
    df = df.loc[valid_ids]
    occ = (df > 2048).astype(int)

    n = occ.shape[0]
    if n == 0:
        print("No contigs found.", file=sys.stderr)
        with open(args.output_file, 'w') as f:
            f.write("contig\tbin\n")
        return

    print(f"Number of contigs: {n}", file=sys.stderr)
    vecs = [occ.iloc[i].values for i in range(n)]

    dist = np.ones((n, n))
    np.fill_diagonal(dist, 0)

    total_pairs = n * (n - 1) // 2
    processed = 0
    for i, j in combinations(range(n), 2):
        a = np.sum((vecs[i] == 1) & (vecs[j] == 1))
        b = np.sum((vecs[i] == 0) & (vecs[j] == 1))
        c = np.sum((vecs[i] == 1) & (vecs[j] == 0))
        d = np.sum((vecs[i] == 0) & (vecs[j] == 0))
        table = [[a, b], [c, d]]
        if a == 0 or b == 0 or c == 0:
            p = 1.0
        else:
            _, p = fisher_exact(table, alternative='two-sided')
        dist[i, j] = p
        dist[j, i] = p
        processed += 1
        if processed % 10000 == 0:
            print(f"Processed {processed}/{total_pairs} pairs", file=sys.stderr)

    if args.tsne:
        print("Computing Spearman correlation and transforming distance matrix (aligning with original code)...", file=sys.stderr)

        # 1. 计算 Spearman 相关系数矩阵 (occ 是 contigs × samples)
        corr_matrix, _ = spearmanr(occ.T)

        # 2. 根据相关系数调整距离矩阵：如果 r(i,j) < 0，则 dist_transformed(i,j) = -dist(i,j)
        dist_transformed = dist.copy()
        neg_mask = corr_matrix < 0
        dist_transformed[neg_mask] = -dist_transformed[neg_mask]

        # 3. 计算 Spearman 距离矩阵
        D = squareform(pdist(dist_transformed, 'correlation'))
        np.fill_diagonal(D, 0)

        print(f"Applying t-SNE to {args.tsne_dim} dimensions with perplexity {args.tsne_perplexity}...", file=sys.stderr)
        tsne = TSNE(
            n_components=args.tsne_dim,
            random_state=2015,
            perplexity=args.tsne_perplexity,
            metric='precomputed',
            init='random'
        )
        X_tsne = tsne.fit_transform(D)  # 使用变换后的距离矩阵 D
        print("t-SNE completed. Clustering with DBSCAN on t-SNE space...", file=sys.stderr)
        clustering = DBSCAN(eps=args.eps, min_samples=args.min_samples)
        labels = clustering.fit_predict(X_tsne)
    else:
        print("Clustering with DBSCAN on original distance matrix...", file=sys.stderr)
        clustering = DBSCAN(eps=args.eps, min_samples=args.min_samples,
                            metric='precomputed')
        labels = clustering.fit_predict(dist)

    contigs = occ.index.tolist()
    bin_map = {}
    for idx, label in enumerate(labels):
        if label == -1:
            bin_map[contigs[idx]] = "unbinned"
        else:
            bin_map[contigs[idx]] = f"bin_{label}"

    cluster_df = pd.DataFrame(bin_map.items(), columns=['contig', 'bin'])
    cluster_df.to_csv(args.output_file, sep='\t', index=False)


    n_bins = len(set(bin_map.values())) - 1
    print(f"Number of bins: {n_bins}", file=sys.stderr)

if __name__ == "__main__":
    main()
