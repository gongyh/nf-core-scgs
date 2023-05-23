#!/usr/bin/env python3

import sys

ctg_genes_file = sys.argv[1]
prokka_tsv = sys.argv[2]

ctg_gene_dict = dict()
with open(ctg_genes_file) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        ctg_gene_dict[cl[1]] = cl[0]

header = True
with open(prokka_tsv) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        if header:
            print("ctg_id\t" + line.strip())
            header = False
            continue
        if cl[1] != "gene":
            ctg_id = ctg_gene_dict[cl[0]]
            print(ctg_id + "\t" + line.strip())
