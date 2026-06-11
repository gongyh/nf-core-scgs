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
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
def main():
    parser = argparse.ArgumentParser(description='Co-occurrence binning')
    parser.add_argument('coverage_file', help='coverage matrix TSV')
    parser.add_argument('output_file', help='output clusters TSV')
    parser.add_argument('--n_clusters', type=int, default=None,
                        help='Number of clusters for single-sample KMeans (default: auto)')
    args = parser.parse_args()

    df = pd.read_csv(args.coverage_file, sep='\t', index_col=0)
    if df.shape[1] != 1:
        print("Error: Input must have exactly one column (single sample).", file=sys.stderr)
        sys.exit(1)

    depths = df.iloc[:, 0].values.reshape(-1, 1)
    scaler = StandardScaler()
    depths_scaled = scaler.fit_transform(depths)

    n_contigs = len(depths_scaled)
    if args.n_clusters is not None:
        n_clusters = args.n_clusters
    else:
        n_clusters = min(20, max(2, n_contigs // 10))

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(depths_scaled)

    contigs = df.index.tolist()
    with open(args.output_file, 'w') as f:
        f.write("contig\tbin\n")
        for contig, label in zip(contigs, labels):
            f.write(f"{contig}\tbin_{label}\n")

    n_bins = len(set(labels))
    print(f"Number of bins: {n_bins}", file=sys.stderr)

if __name__ == "__main__":
    main()
