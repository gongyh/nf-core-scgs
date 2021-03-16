#!/usr/bin/env python3

import sys

if len(sys.argv)!=3:
    print("Usage: python3 kofam_postprocess.py ko_KO.txt test_KOs_mapper.txt")
    exit(0)

ko_KO = sys.argv[1]
mapper = sys.argv[2]

KO_dict = dict()
with open(ko_KO) as fh:
    for line in fh:
        cl = line.strip().split("\t")
        KO = cl[2]
        KO_dict[KO] = line.strip()

header = 1
with open(mapper) as fh:
    for line in fh:
        if header:
            header -= 1
            continue
        cl = line.strip().split("\t")
        if len(cl)==1:
            continue
        KO = cl[1]
        detail = KO_dict[KO].split("\t")
        ko = "ko" + detail[0]
        pathway_detail = ko +" "+ detail[1].split(" [")[0]
        KO_details = detail[3].split("; ")
        KO_detail = detail[2]+" "+KO_details[0]
        KO_definition = KO_details[1]
        print(pathway_detail+"\t"+KO_detail+"\t"+KO_definition)

