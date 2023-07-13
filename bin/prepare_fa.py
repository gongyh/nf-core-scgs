#!/usr/bin/env python
import os
import csv
import shutil

tpath = "../split"
dir_l = os.listdir(tpath)
tsv_list = []
dir_list = []
sample_dict = {}
for cur_file in dir_l:
    cur_path = os.path.join(tpath, cur_file)
    if os.path.isfile(cur_path) and os.path.splitext(cur_file)[1] == ".tsv":
        tsv_list.append(cur_file)
    if os.path.isdir(cur_path) and cur_file.endswith("_Bacteria"):
        dir_list.append(cur_file)

for i in tsv_list:
    sample = i.split("_")[0]
    file_list = []
    fpath = os.path.join(tpath, i)
    with open(fpath, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if float(row[11]) >= 40.0 and float(row[12]) <= 10.0 and row[0] != "no-hit":
                file_list.append(row[0])
    sample_dict[sample] = file_list

if sample_dict:
    for l in dir_list:
        dsample = l.split("_")[0]
        for fname in sample_dict[dsample]:
            base_file = os.path.join(tpath, l, fname + ".fasta")
            target_file = os.path.join("./fa", dsample + "_" + fname + ".fasta")
            shutil.copyfile(base_file, target_file)
