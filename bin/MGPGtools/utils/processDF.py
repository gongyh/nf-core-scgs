import os
import subprocess
import pandas as pd
from MGPGtools.utils.odgi import *
from MGPGtools.utils.common import delete_files
from MGPGtools.utils.gfa import nodeLength


## 处理gene_df的每一行，用于多进程
def processRow(row, genome_df, ref):
    # 取得行
    row = row[1]
    result = {}
    tmp_df = pd.DataFrame(data=row)
    # 每行第3列之后是1的行索引为该基因的node
    column_indexes = tmp_df[3:][tmp_df.iloc[3:] == 1].dropna(axis=0, how="all").index.tolist()
    # 列索引加入"path.name"索引
    column_indexes.insert(0, "path.name")
    # genome dataframe根据基因的列提取部分列
    genome_gene_df = genome_df.loc[:, column_indexes]
    # 遍历path dataframe，获得每个基因组与reference gene相同的node组成的list
    result[row["path.name"]] = {}
    for i, r in genome_gene_df.iterrows():
        if ref in str(r["path.name"]):
            continue
        columns = list(genome_gene_df.columns[2:][r[2:] == 1])
        if len(columns) == 0:
            continue
        else:
            result[row["path.name"]][r["path.name"]] = columns
    return result


def processDf(df, coreGene):
    result = {}
    # df是一个dataframe，前部分行是gene，后部分行是genome paths
    # 提取genome path的行构成新的dataframe
    genome_df = df[df["path.name"].str.contains("GCA|GCF") == True]
    # 提取gene行构成新的dataframe
    gene_df = df[df["path.name"].str.contains("GCA|GCF") == False]
    # 遍历gene dataframe的每一行
    for index, row in gene_df.iterrows():
        # tmp_df = pd.DataFrame(data=row).transpose()
        # 找到从第3列开始(第三列开始是node，前两列分别是path.name, path.conut)所有是1的nodeid，即是该基因的node，获得这些node构成的列索引
        column_indexes = row.columns[3:][row.iloc[0, 3:] == 1].tolist()
        # 列索引加入"path.name"索引
        column_indexes.insert(0, "path.name")
        # genome dataframe根据基因的列提取部分列，不全为0的行是突变区域，全为0的行是非突变区域
        genome_gene_df = genome_df.loc[:, column_indexes]
        genome_gene_df = genome_gene_df[(genome_gene_df.iloc[:, 1:] != 0).any(axis=1)]
        # 遍历基因组path dataframe，获得每个基因组与reference gene不同的node组成的list
        for i, r in genome_gene_df.iterrows():
            columns = list(genome_gene_df.columns[2:][r[2:] == 1])
            result[row["path.name"]][r["path.name"]] = columns
    return result


# def mergeDF(df):
#     df_merged = (
#         df.groupby(df.iloc[:, 0].str.split("#").str[:2].str.join("#"))
#         .apply(merge_rows)
#         .reset_index(drop=True)
#     )
#     genome_df = df_merged[df_merged["path.name"].str.contains("GCA|GCF") == True]
#     # 提取gene行构成新的dataframe
#     gene_df = df_merged[df_merged["path.name"].str.contains("GCA|GCF") == False]
#     print("[genome_df finished****]")
#     print("[gene_df finished****]")
#     return gene_df, genome_df


# def merge_rows(row_group):
#     merged_row = row_group.iloc[0].copy()
#     merged_row.iloc[0] = "#".join(row_group.iloc[0].iloc[0].split("#")[:2])
#     for i in range(1, len(row_group.columns)):
#         merged_row.iloc[i] = row_group.iloc[:, i].sum()
#     return merged_row


# 根据每个基因的位置提取子图
def extractGenesOg(genePath, ogFile, outdir, geneTag, geneLength, genomeListExceptRef):
    geneName = geneTag[genePath]
    gene_length = geneLength[geneName]
    # 每个基因的odgi文件
    extractSortedOg = os.path.join(outdir, genePath + ".sorted.og")
    extractGfa = os.path.join(outdir, genePath + ".gfa")
    extractCmd = [
        "odgi",
        "extract",
        "-i",
        ogFile,
        "-r",
        genePath,
        "-d",
        "3000",
        "-o",
        "-",
    ]
    extractSortCmd = ["odgi", "sort", "-i", "-", "-o", extractSortedOg]
    p1 = subprocess.Popen(extractCmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(extractSortCmd, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p2.communicate()
    ogView(extractSortedOg, extractGfa, 1)
    nodeL = nodeLength(extractGfa)
    if_success, stdout, stderr = ogPath(extractSortedOg, 1)
    lines = stdout.split("\n")
    lines.pop()
    columns = lines[0].split("\t")
    data = [l.split("\t") for l in lines[1:]]
    # 每个基因的矩阵
    df = pd.DataFrame(data, columns=columns)
    cols_to_extract = df.columns[3:][df.iloc[0, 3:] == "1"].tolist()
    cols_to_extract.insert(0, "path.name")
    genomeDF = df[cols_to_extract]
    # genome_df: 记录矩阵中出现的基因组名称
    # variant_genome: 这个基因上突变碱基总数大于0.2的基因组
    genome_df = []
    variant_genome = []
    for index, row in genomeDF.iterrows():
        # 去除参考基因组
        if genePath in str(row["path.name"]):
            continue
        # 基因组名称
        genomeName = row["path.name"].split("#")[0] + "." + row["path.name"].split("#")[1]
        if genomeName not in genome_df:
            genome_df.append(genomeName)
        else:
            if genomeName not in variant_genome:
                continue
        zeroColumns = list(genomeDF.columns[1:][row[1:] == "0"])
        length = 0
        for n in zeroColumns:
            length += nodeL[n[5:]]
        if length / gene_length > 0.2:
            if genomeName not in variant_genome:
                variant_genome.append(genomeName)
    absence_genome = list(set(genomeListExceptRef) - set(genome_df))
    variant_genome.extend(absence_genome)
    absenceGene = {}
    absenceGene[geneName] = variant_genome
    # delete_files(extractSortedOg)
    # delete_files(extractGfa)
    return absenceGene
