#!/usr/bin/env python3

## split SAG assembly according to contig annotation

import os
import sys
from Bio import SeqIO

if len(sys.argv) != 5:
    print("Usage: python splitAss.py RG1.ctg200.fasta RG1.blobDB.bestsum.table.txt family out_dir")
    exit(0)

fa = sys.argv[1]
ann = sys.argv[2]
level = sys.argv[3]
out_dir = sys.argv[4]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

levels = ["superkingdom", "phylum", "order", "family", "genus", "species"]

record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

if level not in levels:
    level = "order"

annCol = -1  # which colum corresponds to the specified level

with open(ann) as fh:
    for line in fh:
        if line[0:2] == "##":  # comment line, skip
            continue
        elif line[0:6] == "# name":  # header line
            cl = line.strip().split("\t")
            for item in cl:
                il = item.split(".")
                if len(il) == 3 and il[0] == level and il[1] == "t":  # eg. order.t.12
                    annCol = int(il[2]) - 1
            continue
        cl = line.strip().split("\t")
        contig_id = cl[0]
        assert annCol > 0
        annotation = cl[annCol]
        fname = annotation.replace(" ", "_") + ".fasta"
        with open(out_dir + "/" + fname, "a") as output_handle:
            SeqIO.write(record_dict[contig_id], output_handle, "fasta")
