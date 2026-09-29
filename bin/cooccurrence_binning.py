#!/usr/bin/env python3
"""
Co-occurrence binning based on Fisher's exact test and DBSCAN.
Input: coverage matrix (TSV: contigs × samples, values = coverage)
Output: clusters.tsv (contig -> bin), coverage and p-value heatmaps,
        and a bin-coloured embedding plot when --tsne is enabled.
Usage: python cooccurrence_binning.py coverage.tsv filtered_ids.txt clusters.tsv
"""

import sys
import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, rankdata
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from itertools import combinations


def plot_heatmap(
    values, rows, columns, output_path, title, color_label, cmap="viridis", vmin=None, vmax=None, cluster=False
):
    fig = plt.figure(figsize=(12, 8), constrained_layout=True)
    if cluster and values.size:
        grid = fig.add_gridspec(2, 3, width_ratios=[1, 6, 0.25], height_ratios=[1, 6])
        ax = fig.add_subplot(grid[1, 1])
        color_ax = fig.add_subplot(grid[1, 2])
        row_order = np.arange(len(rows))
        column_order = np.arange(len(columns))
        for observations, position, orientation in [(values, grid[1, 0], "left"), (values.T, grid[0, 1], "top")]:
            if len(observations) < 2:
                continue
            distances = pdist(observations, metric="euclidean")
            tree = linkage(distances, method="average")
            tree_ax = fig.add_subplot(position)
            if np.any(distances):
                dendrogram(
                    tree,
                    ax=tree_ax,
                    orientation=orientation,
                    no_labels=True,
                    color_threshold=0,
                    above_threshold_color="#444444",
                )
            else:
                tree_ax.text(
                    0.5,
                    0.5,
                    "Identical profiles",
                    ha="center",
                    va="center",
                    rotation=90 if orientation == "left" else 0,
                    transform=tree_ax.transAxes,
                )
            tree_ax.set_axis_off()
            if orientation == "left":
                row_order = leaves_list(tree)
                tree_ax.invert_yaxis()
            else:
                column_order = leaves_list(tree)
        values = values[np.ix_(row_order, column_order)]
        rows = [rows[index] for index in row_order]
        columns = [columns[index] for index in column_order]
    else:
        ax = fig.add_subplot(111)
        color_ax = None
    if cluster and values.size:
        fig.suptitle(title)
    else:
        ax.set_title(title)
    ax.set_xlabel("Contig")
    ax.set_ylabel("Sample" if color_label == "log2(coverage + 1)" else "Contig")
    if values.size == 0:
        ax.text(0.5, 0.5, "No contigs remain after filtering", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        image = ax.imshow(values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
        for axis, names in [(ax.xaxis, columns), (ax.yaxis, rows)]:
            positions = np.unique(np.linspace(0, len(names) - 1, min(len(names), 30)).astype(int))
            axis.set_ticks(positions)
            axis.set_ticklabels([names[i] for i in positions])
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=8)
        plt.setp(ax.get_yticklabels(), fontsize=8)
        fig.colorbar(image, ax=ax, cax=color_ax, label=color_label)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_embedding(embedding, labels, output_path):
    fig = plt.figure(figsize=(10, 8), constrained_layout=True)
    dimensions = embedding.shape[1]
    ax = fig.add_subplot(111, projection="3d" if dimensions >= 3 else None)
    bin_labels = sorted(set(labels))
    colors = plt.get_cmap("tab20")
    for index, label in enumerate(bin_labels):
        points = embedding[labels == label]
        color = "#888888" if label == -1 else colors(index % 20)
        name = "unbinned" if label == -1 else f"bin_{label}"
        coordinates = [points[:, 0], points[:, 1] if dimensions >= 2 else np.zeros(len(points))]
        if dimensions >= 3:
            coordinates.append(points[:, 2])
        ax.scatter(*coordinates, color=color, label=name, s=20, alpha=0.8)
    ax.set_title("Contig t-SNE embedding")
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2" if dimensions >= 2 else "")
    if dimensions >= 3:
        ax.set_zlabel("t-SNE 3")
    if len(bin_labels) <= 20:
        ax.legend(title="Bin", loc="best")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Co-occurrence binning")
    parser.add_argument("coverage_file", help="coverage matrix TSV")
    parser.add_argument("filtered_ids", help="filtered contigs ID list")
    parser.add_argument("output_file", help="output clusters TSV")
    parser.add_argument("--eps", type=float, default=0.05, help="DBSCAN eps (default: 0.05)")
    parser.add_argument("--min_samples", type=int, default=2, help="DBSCAN min_samples (default: 2)")
    parser.add_argument("--tsne", action="store_true", help="Apply t-SNE before DBSCAN")
    parser.add_argument("--tsne_dim", type=int, default=2, help="t-SNE output dimensions (default: 2)")
    parser.add_argument("--tsne_perplexity", type=float, default=30.0, help="t-SNE perplexity (default: 30)")
    args = parser.parse_args()
    if args.tsne_dim < 1:
        parser.error("--tsne_dim must be at least 1")
    if args.tsne_perplexity <= 0:
        parser.error("--tsne_perplexity must be positive")
    output_dir = Path(args.output_file).parent
    contig_ids = []
    with open(args.filtered_ids) as f:
        for line in f:
            contig_id = line.strip()
            if contig_id:
                contig_ids.append(contig_id)
    print(f"Filtered contigs from fasta: {len(contig_ids)}", file=sys.stderr)
    df = pd.read_csv(args.coverage_file, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    valid_ids = [id for id in contig_ids if id in df.index]
    print(f"Valid contigs in matrix: {len(valid_ids)}", file=sys.stderr)
    df = df.loc[valid_ids]
    plot_heatmap(
        np.log2(df.clip(lower=0).to_numpy().T + 1),
        df.columns.tolist(),
        valid_ids,
        output_dir / "coverage_heatmap.png",
        "Filtered contig coverage",
        "log2(coverage + 1)",
        cluster=True,
    )
    if len(valid_ids) == 0:
        print("No valid contigs found after filtering.", file=sys.stderr)
        with open(args.output_file, "w") as f:
            f.write("contig\tbin\n")
        plot_heatmap(
            np.empty((0, 0)),
            [],
            [],
            output_dir / "pvalue_heatmap.png",
            "Contig-contig Fisher exact test p-values",
            "p-value",
            cmap="magma_r",
            vmin=0,
            vmax=1,
        )
        return
    occ = (df > 2048).astype(int)

    n = occ.shape[0]
    if n == 0:
        print("No contigs found.", file=sys.stderr)
        with open(args.output_file, "w") as f:
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
        b = np.sum((vecs[i] == 1) & (vecs[j] == 0))
        c = np.sum((vecs[i] == 0) & (vecs[j] == 1))
        d = np.sum((vecs[i] == 0) & (vecs[j] == 0))
        table = [[a, b], [c, d]]
        _, p = fisher_exact(table, alternative="two-sided")
        dist[i, j] = p
        dist[j, i] = p
        processed += 1
        if processed % 10000 == 0:
            print(f"Processed {processed}/{total_pairs} pairs", file=sys.stderr)

    dist2 = np.maximum(dist, 1e-16)
    plot_heatmap(
        np.log10(dist2),
        valid_ids,
        valid_ids,
        output_dir / "pvalue_heatmap.png",
        "Contig-contig Fisher exact test p-values",
        "p-value (log10)",
        cmap="Greys_r",
        cluster=True,
    )

    if args.tsne and n > 1:
        print(
            "Computing Spearman correlation and transforming distance matrix (aligning with original code)...",
            file=sys.stderr,
        )

        # 1. Compute Spearman correlation matrix (occ is contigs x samples)
        corr_matrix = df.T.corr(method="spearman").to_numpy()

        # 2. Adjust distance matrix based on correlation: if r(i,j) < 0, then dist_transformed(i,j) = -dist(i,j)
        dist_transformed = dist.copy()
        neg_mask = corr_matrix < 0
        dist_transformed[neg_mask] = -dist_transformed[neg_mask]

        # 3. Compute Spearman distance matrix
        dist_ranked = np.apply_along_axis(rankdata, 1, dist_transformed)
        D = squareform(pdist(dist_ranked, "correlation"))
        #D = np.nan_to_num(D, nan=1.0, posinf=1.0, neginf=0.0)
        #D = np.maximum(D, 0)
        np.fill_diagonal(D, 0)

        perplexity = min(args.tsne_perplexity, n - 1)
        print(f"Applying t-SNE to {args.tsne_dim} dimensions with perplexity {perplexity}...", file=sys.stderr)
        tsne = TSNE(
            n_components=args.tsne_dim,
            random_state=2015,
            early_exaggeration=4.0, max_iter=5000,
            perplexity=perplexity,
            metric="precomputed",
            init="random",
            method="barnes_hut" if args.tsne_dim <= 3 else "exact",
        )
        X_tsne = tsne.fit_transform(D)  # use transformed distance matrix D
        print("t-SNE completed. Clustering with DBSCAN on t-SNE space...", file=sys.stderr)
        clustering = DBSCAN(eps=args.eps, min_samples=args.min_samples)
        labels = clustering.fit_predict(X_tsne)
    else:
        print("Clustering with DBSCAN on original distance matrix...", file=sys.stderr)
        clustering = DBSCAN(eps=args.eps, min_samples=args.min_samples, metric="precomputed")
        labels = clustering.fit_predict(dist)

    if args.tsne:
        if n == 1:
            X_tsne = np.zeros((1, args.tsne_dim))
            print("Only one contig; plotting at the origin without fitting t-SNE.", file=sys.stderr)
        plot_embedding(X_tsne, labels, output_dir / "tsne_embedding.png")

    contigs = occ.index.tolist()
    bin_map = {}
    for idx, label in enumerate(labels):
        if label == -1:
            bin_map[contigs[idx]] = "unbinned"
        else:
            bin_map[contigs[idx]] = f"bin_{label}"

    cluster_df = pd.DataFrame(bin_map.items(), columns=["contig", "bin"])
    cluster_df.to_csv(args.output_file, sep="\t", index=False)

    n_bins = len(set(bin_map.values()) - {"unbinned"})
    print(f"Number of bins: {n_bins}", file=sys.stderr)


if __name__ == "__main__":
    main()
