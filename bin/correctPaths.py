#!/usr/bin/env python3
import sys

contigsFasta = sys.argv[1]
contigsPath = sys.argv[2]
contigsID = []

contigsCorrPath = contigsPath.split(".")[0] + ".correct.paths"
with open(contigsFasta, "r") as ffasta:
    fastaLines = ffasta.readlines()
    for fastaLine in fastaLines:
        if fastaLine.startswith(">"):
            ID = fastaLine.strip().replace(">","")
            contigsID.append(ID)


with open(contigsPath, "r") as fcontigs, open(contigsCorrPath, "w") as fcorrect:
    contigsLines = fcontigs.readlines()
    for contigsLine in contigsLines:
        if contigsLine.startswith("NODE"):
            l = contigsLine.strip()
            if l.endswith("'"):
                if l[:-1] in contigsID:
                    fcorrect.write(l+"\n")
                else:
                    for i in contigsID:
                        if i[:i.find("_length")] == l[:l.find("_length")]:
                            fcorrect.write(i+"'\n")
            else:
                if l in contigsID:
                    fcorrect.write(l+"\n")
                else:
                    for i in contigsID:
                        if i[:i.find("_length")] == l[:l.find("_length")]:
                            fcorrect.write(i+"\n")
        else:
            fcorrect.write(contigsLine)
