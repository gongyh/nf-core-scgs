import re
import sys
import os
from MGPGtools.utils.common import *
from MGPGtools.utils.odgi import *
from MGPGtools.utils.meta import get_info, get_chromosome


def gfa_parse_link_path(gfaFile):
    """Parse link and path from gfa files

    Returns:
        gfa_path: paths in gfa files
        gfa_link: links in gfa files
    """
    gfa_path = {}
    gfa_node = {}
    with open(gfaFile, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "S":
                gfa_node["node_" + row[1]] = row[2]
            if row[0] == "P":
                # chrom = re.sub(r"\[.*?\]", "", row[1])
                chrom = row[1]
                gfa_path[chrom] = {}
                if "," in row[2]:
                    gfa_path[chrom]["path"] = row[2].split(",")
                    gfa_path[chrom]["overlap"] = row[3]
                    if int(gfa_path[chrom]["path"][0][:-1]) < int(
                        gfa_path[chrom]["path"][-1][:-1]
                    ):
                        gfa_path[chrom]["tag"] = "forward"
                    else:
                        gfa_path[chrom]["tag"] = "reverse"
                else:
                    gfa_path[chrom]["path"] = [row[2]]
                    gfa_path[chrom]["overlap"] = row[3]
                    if row[2][-1] == "+":
                        gfa_path[chrom]["tag"] = "forward"
                    else:
                        gfa_path[chrom]["tag"] = "reverse"
    return gfa_node, gfa_path


def nodeStr(gfaFile):
    node = {}
    with open(gfaFile, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "S":
                node[row[1]] = row[2]
    return node


def nodeLength(gfaFile):
    nodeLen = {}
    with open(gfaFile, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            r = []
            row = line.strip().split("\t")
            if row[0] == "S":
                nodeLen[row[1]] = len(row[2])
    return nodeLen
