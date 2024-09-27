import re
import sys
import os
from utils.common import *
from collections import deque
from utils.meta import get_info, get_chromosome
from utils.gff import extract_gene


def gfa_parse_link_path(database, name):
    """Parse link and path from gfa files

    Returns:
        gfa_path: paths in gfa files
        gfa_link: links in gfa files
    """
    gfa_path = dict()
    gfa_link = dict()
    unsigned_gfa_link = dict()
    details = get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)
    gfa_file = os.path.join(
        database, "databases", "gfa", details["class"], name + ".gfa"
    )
    with open(gfa_file, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "L":
                prec_node = row[1] + row[2]
                prec_node_number = row[1]
                rear_node = row[3] + row[4]
                rear_node_number = row[3]
                if prec_node in gfa_link:
                    gfa_link[prec_node].append(rear_node)
                    unsigned_gfa_link[prec_node_number].append(rear_node_number)
                else:
                    gfa_link[prec_node] = [rear_node]
                    unsigned_gfa_link[prec_node_number] = [rear_node_number]
            if row[0] == "P":
                # chrom = re.sub(r"\[.*?\]", "", row[1])
                chrom = row[1]
                gfa_path[chrom] = {}
                if "," in row[2]:
                    gfa_path[chrom]["path"] = row[2].split(",")
                    gfa_path[chrom]["overlap"] = row[3]
                else:
                    p = []
                    gfa_path[chrom]["path"] = p.append(row[2])
                    gfa_path[chrom]["overlap"] = row[3]
    return gfa_link, unsigned_gfa_link, gfa_path


def core_genome_extract():
    """extract core sequence/dispensable sequence"""
    ref_nodes = list()
    ref_chrom_list = list()
    ref_chrom_path_list = list()
    alt_path_range = dict()
    info, chrom_size = get_chromosome(database, name)
    details = get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)
    gfa_file = os.path.join(
        database, "databases", "gfa", details["class"], name + ".gfa"
    )
    ref_genome = info["ref"].replace(".", "#")
    # ref_chrom_path_list: path形式list
    # ref_chrom_list: chromosome name形式list
    for chrom in chrom_size.keys():
        ref_chrom_path_list.append(ref_genome + "#" + chrom)
        ref_chrom_list.append(chrom)
    gfa_link, unsigned_gfa_link, gfa_path = gfa_parse_link_path()
    ## k: 染色体
    for k in ref_chrom_path_list:
        ref_nodes.extend(gfa_path[k]["path"])
    unsigned_ref_nodes = [x[:-1] for x in ref_nodes]
    # 除参考基因组外的其他path
    alt_path_range = {
        key: value["path"]
        for key, value in gfa_path.items()
        if key not in ref_chrom_path_list
    }
    unsigned_alt_path_range = {
        key: [x[:-1] for x in value] for key, value in alt_path_range.items()
    }
    nodes_position_list, segment_length = nodes_position(ref_chrom_path_list)
    genes = extract_gene(database, name, ref_chrom_list)
    variant_dict = {}
    variant_nodes_dict = {}
    genes_list = []
    variant_range = []
    if alt_path_range:
        for k, v in alt_path_range.items():
            variant_dict[k] = dict()
            variant_dict[k]["start"] = v[0]
            variant_dict[k]["end"] = v[-1]
            # 初始化variant_nodes_dict
            for node in v:
                if node[:-1] in variant_nodes_dict:
                    variant_nodes_dict[node[:-1]].append(0)
                else:
                    variant_nodes_dict[node[:-1]] = [0]
            # 正向序列
            if v[0].endswith("+"):
                new_index = conversed_nodes_index = [
                    v.index(x) for x in v if x in ref_nodes
                ]
            else:
                rev_ref_nodes = []
                for n in ref_nodes:
                    new_element = n.replace("+", "-")
                    rev_ref_nodes.append(new_element)
                new_index = conversed_nodes_index = [
                    v.index(x) for x in v if x in rev_ref_nodes
                ]
            for i in range(len(conversed_nodes_index) - 1):
                if conversed_nodes_index[i + 1] - conversed_nodes_index[i] > 1:
                    var_bases = 0
                    start_index = conversed_nodes_index[i]
                    end_index = conversed_nodes_index[i + 1]
                    sub_list = v[(start_index + 1) : end_index]
                    for j in sub_list:
                        var_bases += segment_length[j[:-1]]
                    if var_bases <= 3:
                        fill_numbers = list(
                            range(
                                conversed_nodes_index[i] + 1,
                                conversed_nodes_index[i + 1],
                            )
                        )
                        new_index.extend(fill_numbers)
            nodes_list = [v[i][:-1] for i in new_index]
            for n in nodes_list:
                variant_nodes_dict[n][-1] = 1
        common_nodes1 = [
            key
            for key, value in variant_nodes_dict.items()
            if all(val == 1 for val in value)
        ]

        check_directory(os.path.join(outdir, "tmp"))
        # 变异区域范围与ref范围对应
        position_file_path = os.path.join(outdir, "tmp", "position.txt")
        write_positions_file(variant_dict, position_file_path)
        og_file_path = os.path.join(outdir, "tmp", name + ".og")
        if len(ref_chrom_list) == 1:
            ODGIBuildCmd = (
                ["odgi", "build", "-g", gfa_file, "-o", og_file_path]
                if threads == 1
                else [
                    "odgi",
                    "build",
                    "-g",
                    gfa_file,
                    "-o",
                    og_file_path,
                    "-t",
                    threads,
                ]
            )
            is_success, stdout, stderr = run(ODGIBuildCmd)
            ODGIPositionCmd = (
                [
                    "odgi",
                    "position",
                    "-i",
                    og_file_path,
                    "-G",
                    position_file_path,
                    "-r",
                    ref_chrom_path_list[0],
                ]
                if threads == 1
                else [
                    "odgi",
                    "position",
                    "-i",
                    og_file_path,
                    "-G",
                    position_file_path,
                    "-r",
                    ref_chrom_path_list[0],
                    "-t",
                    threads,
                ]
            )
            is_success, stdout, stderr = run(ODGIPositionCmd)
            output = stdout.strip().split("\n")
            for output_line in output:
                if output_line[0] == "#":
                    continue
                node = output_line.strip().split("\t")[0].split(",")[0]
                target_p = output_line.strip().split("\t")[1].split(",")[1]
                dist = output_line.strip().split("\t")[2]
                strand = output_line.strip().split("\t")[3]
                for v in variant_dict.values():
                    if v["start"][:-1] == node:
                        if strand == "+":
                            v["start_position"] = int(target_p) + int(dist)
                        else:
                            v["start_position"] = int(target_p) - int(dist)
                    if v["end"][:-1] == node:
                        if strand == "+":
                            v["end_position"] = int(target_p) + int(dist)
                        else:
                            v["end_position"] = int(target_p) - int(dist)
            for v in variant_dict.values():
                if v["start_position"] > v["end_position"]:
                    tmp = v["start_position"]
                    v["start_position"] = v["end_position"]
                    v["end_position"] = tmp
                variant_range.append([v["start_position"], v["end_position"]])
            variant_range = get_merged_ranges(variant_range)
            variant_range_index = []
            for r in variant_range:
                for i in nodes_position_list:
                    if i[0] == i[1]:
                        if r[0] == i[0]:
                            start_node_index = nodes_position_list.index(i)
                        if r[1] == i[0]:
                            end_node_index = nodes_position_list.index(i)
                    else:
                        if i[0] <= r[0] < i[1]:
                            start_node_index = nodes_position_list.index(i)
                        if i[0] < r[1] <= i[1]:
                            end_node_index = nodes_position_list.index(i)
                        if r[0] == i[1]:
                            start_node_index = nodes_position_list.index(i) + 1
                        if r[1] == i[0]:
                            end_node_index = nodes_position_list.index(i) - 1
                variant_range_index.append([start_node_index, end_node_index])
            common_nodes2 = remove_elements(unsigned_ref_nodes, variant_range_index)
            common_nodes = common_nodes1 + common_nodes2
            core_gfa_file_path = os.path.join(outdir, "tmp", name + "_core.gfa")
            set1 = set(variant_nodes_dict.keys())
            set2 = set(common_nodes2)
            if set1.intersection(set2):
                l = set1.intersection(set2)
                print(l)
            else:
                print("No Intersection!")
            core_gfa_file(common_nodes, gfa_file, core_gfa_file_path)


def nodes_position(chrom_path_list):
    s_len = dict()
    range_list = list()
    info, chrom_size = get_chromosome(database, name)
    gfa_link, unsigned_gfa_link, gfa_path = gfa_parse_link_path()
    for k in chrom_path_list:
        range_list.extend(gfa_path[k]["path"])
    gfa_path = os.path.join(
        database,
        "databases",
        "gfa",
        get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)["class"],
        name + ".gfa",
    )
    with open(gfa_path, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            r = list()
            row = line.strip().split("\t")
            if row[0] == "S":
                s_len[row[1]] = len(row[2])
            if row[0] == "L" or row[0] == "P":
                break
    a = []
    start = 1
    for i in range_list:
        if start == 1:
            start = s_len[i[:-1]]
            r = [1, start]
        else:
            r = [start + 1]
            start = start + s_len[i[:-1]]
            r.append(start)
        a.append(r)
    return a, s_len


def get_merged_ranges(ranges):
    if not ranges:
        return []

    sorted_ranges = sorted(ranges)
    merged_ranges = [sorted_ranges[0]]

    for cur_range in sorted_ranges[1:]:
        if merged_ranges[-1][1] >= cur_range[0]:
            merged_ranges[-1][1] = max(merged_ranges[-1][1], cur_range[1])
        else:
            merged_ranges.append(cur_range)

    concat_ranges = []
    for i in range(len(merged_ranges)):
        if i == 0 or merged_ranges[i][0] != merged_ranges[i - 1][1] + 1:
            concat_ranges.append(merged_ranges[i])
        else:
            concat_ranges[-1][1] = merged_ranges[i][1]
    return concat_ranges


def write_positions_file(variant_dict, tmp_path):
    with open(tmp_path, "w") as f:
        for v in variant_dict.values():
            if v["start"].endswith("+"):
                f.write(v["start"] + ",0,+\n")
            if v["start"].endswith("-"):
                f.write(v["start"] + ",0,-\n")
            if v["end"].endswith("+"):
                f.write(v["end"] + ",0,-\n")
            if v["end"].endswith("-"):
                f.write(v["end"] + ",0,+\n")


def combine_chrom_gfa_file(
    ref, chrom_list, chrom_path_list, gfa_path, combine_gfa_path
):
    ref_path_segment = dict()
    ref_segment = ""
    with open(gfa_path, "r") as fin, open(combine_gfa_path, "w") as fout:
        lines = fin.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            row = line.strip().split("\t")
            next_line = lines[i + 1].strip()
            if line[0] == "P" and next_line[0] == "L":
                for j in chrom_path_list:
                    ref_segment += ref_path_segment[j]
                fout.write(line + "\n")
                fout.write("P\t" + ref + "\t" + ref_segment + "\t" + "*\n")
            if row[0] == "P" and ref in row[1]:
                if "[" in row[1]:
                    ref_path_key = re.sub("\[.*?\]", "", row[1])
                    ref_path_segment[ref_path_key] = row[2]
                else:
                    ref_path_segment[row[1]] = row[2]
            else:
                fout.write(line + "\n")


def core_gfa_file(common_nodes, gfa_path, core_gfa_path):
    with open(gfa_path, "r") as fin, open(core_gfa_path, "w") as fout:
        lines = fin.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            row = line.strip().split("\t")
            if row[0] == "H":
                fout.write(line + "\n")
            if row[0] == "S":
                if row[1] in common_nodes:
                    fout.write(line + "\n")
            if row[0] == "L":
                if row[1] in common_nodes or row[3] in common_nodes:
                    fout.write(line + "\n")


def remove_elements(lst1, lst2):
    empty_list = []
    for sublist in lst2:
        start_index, end_index = sublist
        empty_list.extend(lst1[start_index : end_index + 1])
    new_list = [element for element in lst1 if element not in empty_list]
    return new_list
