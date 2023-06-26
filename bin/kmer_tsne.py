#!/usr/bin/env python3

import sys
import h5py
from openTSNE import TSNE, TSNEEmbedding, affinity, initialization

# from openTSNE.callbacks import ErrorLogger
import numpy as np
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: python kmer_tsne.py 4mer.hdf5 out_tsne.tsv [8]")
    exit(0)

kmer_hdf = sys.argv[1]  # from kPAL count
outfile = sys.argv[2]  # e.g. tsne.tsv

cores = 8
if len(sys.argv) >= 4:
    cores = sys.argv[3]  # cpu threads

## 1. first read h5 and convert to matrix
f = h5py.File(kmer_hdf, "r")
fp = f["profiles"]
ctgs = fp.keys()
total = len(ctgs)
profile = np.ones((total, 256), dtype=np.float64)
i = 0
for ctg in ctgs:
    raw = fp[ctg][:]
    norm = raw * 256.0 / sum(raw)
    profile[i,] = norm
    i = i + 1

## 2. Run t-SNE
seed = 211
threads = int(cores)

affinities_multiscale_mixture = affinity.Multiscale(
    profile,
    perplexities=[30, 300],
    metric="cosine",
    n_jobs=threads,
    random_state=seed,
)

init = initialization.pca(profile, random_state=seed)

embedding = TSNEEmbedding(
    init,
    affinities_multiscale_mixture,
    negative_gradient_method="fft",
    n_jobs=threads,
)

embedding1 = embedding.optimize(n_iter=250, exaggeration=6, momentum=0.5)
embedding2 = embedding1.optimize(n_iter=750, exaggeration=1, momentum=0.8)
embedding_multiscale = embedding2.view(np.ndarray)

## 3. Save matrix
df = pd.DataFrame(embedding_multiscale, columns=["Dim1", "Dim2"], index=ctgs)
df.to_csv(outfile, sep="\t", float_format="%.3f")
