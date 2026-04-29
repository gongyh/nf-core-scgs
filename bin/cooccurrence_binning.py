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
from sklearn.cluster import DBSCAN
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(description='Co-occurrence binning')
    parser.add_argument('coverage_file', help='coverage matrix TSV')
    parser.add_argument('output_file', help='output clusters TSV')
    parser.add_argument('--eps', type=float, default=0.05,
                        help='DBSCAN eps (default: 0.05)')
    parser.add_argument('--min_samples', type=int, default=2,
                        help='DBSCAN min_samples (default: 2)')
    args = parser.parse_args()

    df = pd.read_csv(args.coverage_file, sep='\t', index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    occ = (df > 0).astype(int)

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

    print("Clustering with DBSCAN...", file=sys.stderr)
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
