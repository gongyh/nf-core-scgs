#!/usr/bin/env python

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import numpy as np
import pandas as pd
import yaml
import glob
import shutil

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

# Define the path pattern for sample files
file_paths = glob.glob("tda/*.TDA_genus.txt")

# List to store individual sample DataFrames
sample_dfs = []

# Dict to store the most abundant genus of each sample
sample_genus = {}

for fp in file_paths:
    # Read the txt file
    df = pd.read_csv(fp, sep="\\t")

    # Extract sample ID from the file name
    sample_id = fp.split("/")[-1].split(".TDA_genus")[0]

    # Add a column to track the sample ID
    df["sample_id"] = sample_id

    # Append to the list
    sample_dfs.append(df)

    # Find the genus with the highest abundance
    top_genus = df.loc[df["abundance"].idxmax(), "genus"]

    # Step 4: Add to the dictionary
    sample_genus[sample_id] = top_genus

# Combine all samples into one DataFrame (long format)
combined_long = pd.concat(sample_dfs, ignore_index=True)

# Pivot: rows = samples, columns = genera, values = abundance (fill missing with 0)
merged_table = combined_long.pivot_table(
    index="sample_id",  # Rows: sample IDs
    columns="genus",  # Columns: genus names
    values="abundance",  # Values: abundance
    aggfunc="sum",  # Sum abundances if a genus is duplicated in the same sample
    fill_value=0,  # Fill missing genera in a sample with 0
)

adata = sc.AnnData(
    X=merged_table.values,  # Abundance matrix (samples × genera)
    obs=pd.DataFrame(index=merged_table.index),  # Sample metadata (index = sample IDs)
    var=pd.DataFrame(index=merged_table.columns)  # Genus metadata (index = genus names)
)

# Add sample metadata
df_top_genus = pd.DataFrame(  # Convert dict to DataFrame
    sample_genus.values(),  # Convert dict items to list of (key, value) tuples
    columns=["sample_genus"],  # Define column names
    index=sample_genus.keys()
)

adata.obs = adata.obs.join(df_top_genus)  # Merge with existing sample metadata

# Perform clustering
# Build the neighborhood graph (required for both Louvain and Leiden)
n_samples = adata.n_obs
sc.tl.pca(adata, n_comps=min(n_samples // 2, 50)) # Compute PCA (n_comps ≥ n_pcs)
sc.pp.neighbors(
    adata,
    n_neighbors=int(np.sqrt(n_samples)),  # Number of nearest neighbors
    n_pcs=min(n_samples // 2, 50),  # Use PCA components
    use_rep="X_pca"  # Use PCA-reduced data (default)
)
# Run Leiden clustering
sc.tl.leiden(
    adata,
    resolution=1.0,  # Same granularity control as Louvain
    key_added="leiden"  # Store results in adata.obs["leiden"]
)

# Perform umap
sc.tl.umap(adata, random_state=0)
# Round to 10 decimal places to ensure stable hashes
adata.obsm["X_umap"] = np.round(adata.obsm["X_umap"], 10)

# Save results
adata.write_h5ad(f"umap.h5ad")
df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
df.to_pickle(f"umap.pkl")

# Plot clusters
sc.pl.umap(adata, color=["sample_genus", "leiden"], ncols=2, title=["UMAP (Genus)", "UMAP (Leiden)"], save=".pdf")
shutil.copy("figures/umap.pdf", "umap.pdf")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
