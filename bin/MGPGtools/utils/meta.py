import os
import csv
from MGPGtools.utils.taxonomy import Taxonomy
from MGPGtools.utils.common import run


def rank_list(meta, rank_name):
    rankList = []
    with open(meta, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for i in reader:
            index = Taxonomy.rank_index[rank_name]
            rankList.append(i[index])
    rList = list(set(rankList))
    return rList


def get_info(meta, name):
    """
    从Metadata.tsv中获取genus的信息
    """
    k = {
        "domain": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "nodes": "",
        "links": "",
        "genomesNum": "",
        "ref": "",
    }
    with open(meta, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for i in reader:
            if i[5] == name:
                results = [i[0], i[1], i[2], i[3], i[4], i[6], i[7], i[8], i[9]]
                info = dict(zip(k.keys(), results))
                break
    return info


def get_chromosome(database, name):
    chromSize = dict()
    info = get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)
    genomePath = os.path.join(
        database,
        "databases",
        "ref_genome",
        info["class"],
        info["ref"] + "_genomic.fna.gz",
    )
    cmd = ["seqkit", "fx2tab", "-j", "20", "-l", "-n", "-i", "-H", genomePath]
    if_success, stdout, stderr = run(cmd)
    a = stdout.strip().split("\n")
    del a[0]
    for i in a:
        r = i.strip().split("\t")
        chromSize[r[0]] = r[1]
    return info, chromSize


def getPanTxt(database, name):
    genomeList = []
    className = get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)[
        "class"
    ]
    panTxt = os.path.join(
        database,
        "databases",
        "panTxt",
        className,
        name + ".txt",
    )
    with open(panTxt, "r") as file:
        lines = file.readlines()
    for line in lines:
        columns = line.strip().split("\t")
        genomeList.append(columns[0])
    return genomeList
