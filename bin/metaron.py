#!/usr/bin/env python
"""
MetaRon - Metagenomic opeRon Prediction pipeline
Created on Tuesday July 12 14:32:15 2018
Version 1.0
@author: Syed Shujaat Ali Zaidi
"""

from itertools import groupby
from optparse import OptionParser
from inspect import getsourcefile
from os.path import abspath
import csv, re, sys, socket, string, random, os, subprocess, time, shutil, errno, os.path, argparse

start_time = time.time()


def main():
    global process_selector, sample_name, scaf_file, gene_file, paired_merged, output_dir, out_dir, read_type, read_length, raw_read1, raw_read2
    """
    This is metagenomic operon prediction software
    """

    working_dir = os.getcwd()
    ##FOLLOWING PARAMETERS ARE REQUIRED FOR PROGRAM
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--sample", help="Sample name without any dot, underscore or dash")
    parser.add_argument(
        "-p",
        "--process",
        help="1. ago: assembly gene prediction and operon prediciton 2. op: operon prediction only. If 'ago', please provide the following parameters: -n,-rl,-rt,[-pe1,pe2|-pm],",
    )
    parser.add_argument(
        "-rt",
        "--read_type",
        help="Enter read type. 'merge' if the reads are paired-end in two file. 'paired' if the reads are paired-end in one file.",
    )
    parser.add_argument(
        "-rl",
        "--read_length",
        help="Enter 'l'if read length is longer than 128 bases and 'r' if read length is shorter than 128 bases",
    )
    parser.add_argument("-pe1", "--paired_1", help="Enter enter paired read file 1")
    parser.add_argument("-pe2", "--paired_2", help="Enter enter paired read file 2")
    parser.add_argument(
        "-pm", "--paired_merged", help="Enter the paired end read file if both pairedend reads are in one file"
    )
    parser.add_argument(
        "-i", "--igp", help="Select the gene prediction .tab file generated via MetageneMark or Prodigal"
    )
    parser.add_argument("-j", "--isc", help="Select the file containing all scaftigs")
    parser.add_argument("-t", "--tool", help="Enter 1 for MetaGeneMark, 2 for Prodigal")
    parser.add_argument("-o", "--output", help="Enter output destination folder")
    args = parser.parse_args()

    ##ASSIGNING VALUES
    if args.sample:
        sample_name = args.sample
        if type(sample_name) == str:
            ##            print 'OK1'
            pass
        else:
            print("Enter a valid sample name without dot, underscore or dash")
            sys.exit(1)
    if args.process:
        process_selector = args.process
        if process_selector == "ago" or process_selector == "op":
            ##            print 'OK2'
            pass
        else:
            print(
                " Argument_option1: 'ago' for raw data assembly, gene prediction and operon prediciton. Argument_option2: 'op' for operon prediction only"
            )
            sys.exit(1)

    ## if process is "ago"then specific parameters will be available
    if process_selector == "ago":

        if args.read_length:
            read_length = args.read_length
            if type(read_length) == str and read_length == "l" or read_length == "r":
                ##                print 'OK3'
                pass
            else:
                print(
                    "Enter 'l' if read length is longer than 128 bases and 'r' if read length is shorter than 128 bases"
                )
                sys.exit(1)

        ## selection of parameters based on read type "merge"or "paired"
        if args.read_type:
            read_type = args.read_type
            if type(read_type) == str and read_type == "merge" or read_type == "paired":
                ##                print 'OK4'
                pass
                if read_type == "merge":
                    if args.paired_1:
                        raw_read1 = args.paired_1
                        if type(raw_read1) == str and os.path.exists(raw_read1) == True:
                            ##                            print 'OK5'
                            pass
                        else:
                            print(" paired end raw read file 1 path is not correct")
                            sys.exit(1)
                    if args.paired_2:
                        raw_read2 = args.paired_2
                        if type(raw_read2) == str and os.path.exists(raw_read2) == True:
                            ##                            print 'OK6'
                            pass
                        else:
                            print(" paired end raw read file 2 path is not correct")
                            sys.exit(1)
                elif read_type == "paired":
                    if args.paired_merged:
                        paired_merged = args.paired_merged
                        if type(paired_merged) == str and os.path.exists(paired_merged) == True:
                            ##                            print 'OK5'
                            pass
                        else:
                            print(" paired end raw read file path is not correct")
                            sys.exit(1)
            else:
                print(
                    "Read type not correct. 'merge' for paired end reads in two seperate files. 'paired' for paired end reads in one file."
                )
                sys.exit(1)
        if args.output:
            out_dir = args.output
            if type(out_dir) == str and os.path.exists(out_dir) == True:
                print("All parameters checked")
                gene_pred_tool = args.process
                gene_pred_tool = "1"
                pass
            else:
                print("output path is not correct")
                sys.exit(1)

    ## if process is "op"then specific parameters will be available
    elif process_selector == "op":
        if args.igp:
            gene_file = args.igp
            gene_pred_file = args.igp
            if type(gene_pred_file) == str and os.path.exists(gene_pred_file) == True:
                ##                print 'OK1'
                pass
            else:
                print("gene prediction file path is not correct")
                sys.exit(1)
        if args.isc:
            scaf_file = args.isc
            if type(scaf_file) == str and os.path.exists(scaf_file) == True:
                old_pathname = os.path.abspath(args.isc)
                ##                print 'pathname is:  ', old_pathname
                pass
            else:
                print("Scaftig file path is not correct")
                sys.exit(1)
        if args.tool:
            gene_pred_tool = args.tool
            if gene_pred_tool == "1" or gene_pred_tool == "2":
                ##                print 'OK3'
                pass
            else:
                print("1 for MetaGeneMark, 2 for Prodigal")
                sys.exit(1)
        if args.output:
            out_dir = args.output
            if type(out_dir) == str and os.path.exists(out_dir) == True:
                ##                print 'OK4'

                pass
            else:
                print("output path is not correct")
                sys.exit(1)
            print("All parameters checked")

    else:
        print("Please select the right parameter for process")
        sys.exit(1)
    folder_title = sample_name
    print("Sample name:    ", folder_title)

    ##CHECKING IF THE FOLDER ALREADY EXIST THEN REMOVE AND CREATE AGAIN
    if os.path.exists(out_dir):
        tmp = os.path.exists(out_dir + "/MetaRon_" + folder_title)
        if tmp == True:
            pass
        else:
            os.makedirs(out_dir + "/MetaRon_" + folder_title)
        output_dir = out_dir + "MetaRon_" + folder_title
        print("OUTPUT DIRECTORY:    ", output_dir)
    elif not os.path.exists(out_dir):
        tmp = "os.path.exists(out_dir+'MetaRon_'+folder_title)"
        if tmp == True:
            pass
        else:
            shutil.rmtree(out_dir + "MetaRon_" + folder_title)
            os.makedirs(out_dir + "MetaRon_" + folder_title)
        output_dir = out_dir + "MetaRon_" + folder_title
    config_file_check()

    ## REMOVING TEMPS....
    try:
        os.remove(output_dir + "/Gene_coord.tab")
        os.remove(output_dir + "/DIR_Proximons_hdr.tab")
        os.remove(output_dir + "/DIR_Proximons_temp.tab")
        os.remove(output_dir + "/Downstream_seq.fasta")
        os.remove(output_dir + "/Upstream_seq.txt")
        os.remove(output_dir + "/Gene_file_edited.tab")
        os.remove(output_dir + "/Gene_file_detailed.tab")
        os.remove(output_dir + "/Gene_file_UDlen.tab")
        os.remove(output_dir + "/promoter_prediction.txt")
        os.remove(output_dir + "/Promoter_1.tab")
        os.remove(output_dir + "/Promoter_edited.tab")
        os.remove(output_dir + "/IGD_PROM_CLUSTER_noHDR.tab")
        os.remove(output_dir + "/IGD_PROM_RES.tab")
        os.remove(output_dir + "/Operon_File.tab")
        os.remove(output_dir + "/PROM_CLUSTER_noHDR.tab")
        os.remove(output_dir + "/Proximons_hdr.tab")
    except OSError:
        pass

    if process_selector == "ago":
        print(
            "********************************** WHOLE METAGENOME PROCESSING AND OPERON PREDICTION **********************************"
        )
        assembly_process()
        header_replacer_F()
        gene_prediction(scaftig_file, sample_name)
        gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict = data_extraction(
            gene_file, gene_pred_tool
        )
        gene_name, strand, gene_st, gene_end, scaf_name_len_dict = seq_info2(
            gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict
        )
        (
            upstream_seq_len,
            scaftig_name_new,
            gene_name_new,
            strand_new,
            gene_st_new,
            gene_end_new,
            scaf_name4dict_new,
            upstream_st_pos,
            upstream_en_pos,
        ) = upstream_coordinates_extraction(scaftig_name, gene_name, strand, gene_st, gene_end, scaf_name4dict)
        leng_scaftig = scaftig_length(scaf_name4dict, scaf_name_len_dict)
        downstream_seq_len, downstream_st_pos, downstream_en_pos = downstream_coordinate_region(
            scaftig_name, gene_st, gene_end, leng_scaftig, scaftig_name_new, gene_end_new
        )
        saving_data(
            scaftig_name_new, gene_name_new, strand_new, gene_st_new, gene_end_new, upstream_seq_len, downstream_seq_len
        )
        saving_data2(
            scaftig_name_new,
            gene_name_new,
            strand_new,
            gene_st_new,
            gene_end_new,
            upstream_seq_len,
            downstream_seq_len,
            upstream_st_pos,
            upstream_en_pos,
            downstream_st_pos,
            downstream_en_pos,
        )
        (
            newDSS_st_pos,
            newDSS_en_pos,
            newUPS_st_pos,
            newUPS_en_pos,
            gene_name4cord_DSS,
            scaf_name4cord_DSS,
            scaf_name4cord_UPS,
            gene_name4cord_UPS,
        ) = UPS_DSS_Slicing(
            downstream_seq_len,
            downstream_st_pos,
            downstream_en_pos,
            scaftig_name_new,
            gene_name_new,
            upstream_seq_len,
            upstream_st_pos,
            upstream_en_pos,
        )
        downstream_coordinate_file, upstream_coordinate_file = Data_saving3(
            newDSS_st_pos,
            newDSS_en_pos,
            newUPS_st_pos,
            newUPS_en_pos,
            gene_name4cord_DSS,
            scaf_name4cord_DSS,
            scaf_name4cord_UPS,
            gene_name4cord_UPS,
            scaftig_name_new,
            gene_name_new,
            strand_new,
            gene_st_new,
            gene_end_new,
        )
        soucedic = getsouce()
        upstream_seq_file = getgenstring_ups(soucedic, upstream_coordinate_file)
        downstream_seq_file = getgenstring_dss(soucedic, downstream_coordinate_file)
        (
            new_k,
            new_scaftig_no,
            new_gene1_name,
            new_gene2_name,
            new_strand_1,
            new_strand_2,
            IGD,
            IGDclstr,
            new_gene1_clstr,
            new_gene2_clstr,
            new_strand_1clstr,
            new_strand_2clstr,
            new_scaftig_no_clstr,
            new_kclstr,
        ) = IGD_calc(scaftig_name_new, strand_new, gene_name_new, gene_st_new, gene_end_new)
        Data_saving4(
            new_k,
            new_scaftig_no,
            new_gene1_name,
            new_gene2_name,
            new_strand_1,
            new_strand_2,
            IGD,
            IGDclstr,
            new_gene1_clstr,
            new_gene2_clstr,
            new_strand_1clstr,
            new_strand_2clstr,
            new_scaftig_no_clstr,
            new_kclstr,
        )
        promoter_prediction(upstream_seq_file)
        scf_nm, gene_nm, prom_res, counter = Promoter_file_parse()
        Prom_IGD_Clustering()
        Prom_clustering()
        Processing_details(file_name, start_time, counter)

    elif process_selector == "op":
        print(
            "\n\n********************************** METAGENOMIC OPERON PREDICTION **********************************\n\n"
        )
        scaf_len(old_pathname)
    ##        header_replacer_F()
    gff2tab(sample_name)
    gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict = data_extraction(
        gene_file, gene_pred_tool
    )
    gene_name, strand, gene_st, gene_end, scaf_name_len_dict = seq_info2(
        gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict
    )
    (
        upstream_seq_len,
        scaftig_name_new,
        gene_name_new,
        strand_new,
        gene_st_new,
        gene_end_new,
        scaf_name4dict_new,
        upstream_st_pos,
        upstream_en_pos,
    ) = upstream_coordinates_extraction(scaftig_name, gene_name, strand, gene_st, gene_end, scaf_name4dict)
    leng_scaftig = scaftig_length(scaf_name4dict, scaf_name_len_dict)
    ####        downstream_seq_len, downstream_st_pos, downstream_en_pos = downstream_coordinate_region(scaftig_name, gene_st, gene_end, leng_scaftig, scaftig_name_new, gene_end_new)
    downstream_seq_len, downstream_st_pos, downstream_en_pos = downstream_coordinate_region(
        scaftig_name, gene_st, gene_end, leng_scaftig, scaftig_name_new, gene_end_new
    )
    saving_data(
        scaftig_name_new, gene_name_new, strand_new, gene_st_new, gene_end_new, upstream_seq_len, downstream_seq_len
    )
    saving_data2(
        scaftig_name_new,
        gene_name_new,
        strand_new,
        gene_st_new,
        gene_end_new,
        upstream_seq_len,
        downstream_seq_len,
        upstream_st_pos,
        upstream_en_pos,
        downstream_st_pos,
        downstream_en_pos,
    )
    (
        newDSS_st_pos,
        newDSS_en_pos,
        newUPS_st_pos,
        newUPS_en_pos,
        gene_name4cord_DSS,
        scaf_name4cord_DSS,
        scaf_name4cord_UPS,
        gene_name4cord_UPS,
    ) = UPS_DSS_Slicing(
        downstream_seq_len,
        downstream_st_pos,
        downstream_en_pos,
        scaftig_name_new,
        gene_name_new,
        upstream_seq_len,
        upstream_st_pos,
        upstream_en_pos,
    )
    downstream_coordinate_file, upstream_coordinate_file = Data_saving3(
        newDSS_st_pos,
        newDSS_en_pos,
        newUPS_st_pos,
        newUPS_en_pos,
        gene_name4cord_DSS,
        scaf_name4cord_DSS,
        scaf_name4cord_UPS,
        gene_name4cord_UPS,
        scaftig_name_new,
        gene_name_new,
        strand_new,
        gene_st_new,
        gene_end_new,
    )
    soucedic = getsouce_op(pathname)
    upstream_seq_file = getgenstring_ups(soucedic, upstream_coordinate_file)
    downstream_seq_file = getgenstring_dss(soucedic, downstream_coordinate_file)
    (
        new_k,
        new_scaftig_no,
        new_gene1_name,
        new_gene2_name,
        new_strand_1,
        new_strand_2,
        IGD,
        IGDclstr,
        new_gene1_clstr,
        new_gene2_clstr,
        new_strand_1clstr,
        new_strand_2clstr,
        new_scaftig_no_clstr,
        new_kclstr,
    ) = IGD_calc(scaftig_name_new, strand_new, gene_name_new, gene_st_new, gene_end_new)
    Data_saving4(
        new_k,
        new_scaftig_no,
        new_gene1_name,
        new_gene2_name,
        new_strand_1,
        new_strand_2,
        IGD,
        IGDclstr,
        new_gene1_clstr,
        new_gene2_clstr,
        new_strand_1clstr,
        new_strand_2clstr,
        new_scaftig_no_clstr,
        new_kclstr,
    )
    promoter_prediction(upstream_seq_file)
    scf_nm, gene_nm, prom_res, counter = Promoter_file_parse()
    Prom_IGD_Clustering()
    Prom_clustering()
    Processing_details(file_name, start_time, counter)


def assembly_process():
    """The raw reads will be assembled via IDBA."""

    global scaftig_file_unprocessed, scaftig_length_file, sample_name
    print("Assembly by IDBA")

    cmd1 = "cd " + output_dir
    if read_type == "merge":
        cmd2 = "fq2fa --" + read_type + " " + raw_read1 + " " + raw_read2 + " " + output_dir + "/" + sample_name + ".fa"
        if read_length == "l":
            cmd3 = (
                "idba --pre_correction -l "
                + output_dir
                + "/"
                + sample_name
                + ".fa"
                + " -o "
                + output_dir
                + "/"
                + sample_name
                + ".assembly"
            )
        else:
            cmd3 = (
                "idba --pre_correction -r "
                + output_dir
                + "/"
                + sample_name
                + ".fa"
                + " -o "
                + output_dir
                + "/"
                + sample_name
                + ".assembly"
            )
    elif read_type == "paired":
        cmd2 = "fq2fa --" + read_type + " " + paired_merged + " " + output_dir + "/" + sample_name + ".fa"
        if read_length == "l":
            cmd3 = (
                "idba --pre_correction -l "
                + output_dir
                + "/"
                + sample_name
                + ".fa"
                + " -o "
                + output_dir
                + "/"
                + sample_name
                + ".assembly"
            )
        else:
            cmd3 = (
                "idba --pre_correction -r "
                + output_dir
                + "/"
                + sample_name
                + ".fa"
                + " -o "
                + output_dir
                + "/"
                + sample_name
                + ".assembly"
            )
    bs = "/"
    cmd4 = (
        "cp "
        + output_dir
        + "/"
        + sample_name
        + ".assembly/scaffold.fa "
        + output_dir
        + "/"
        + sample_name
        + ".scaffold.fa"
    )  # copy and rename scaffold file to output directory
    cmd5 = (
        "cat "
        + output_dir
        + bs
        + sample_name
        + '.scaffold.fa | awk \'$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\\tlength="; } $0 !~ ">" {c+=length($0);} END { print c; }\' > '
        + output_dir
        + bs
        + sample_name
        + "_len.txt && sed -i '/^$/d' "
        + output_dir
        + bs
        + sample_name
        + "_len.txt"
    )

    ## calculating length of each scaffold

    ##    cmds = [cmd4,cmd5]
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5]
    for cmd in cmds:
        subprocess.call(cmd, shell=True)
    scaftig_file_unprocessed = output_dir + "/" + sample_name + ".scaffold.fa"
    scaftig_length_file = output_dir + "/" + sample_name + "_len.txt"
    print("assembly completed")
    return (scaftig_file_unprocessed, scaftig_length_file, sample_name)


def scaf_len(old_pathname):
    print("Formatting assembly file")
    global pathname
    bs = "/"
    cmd1 = (
        "cat "
        + old_pathname
        + ' | awk \'$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\\tlength="; } $0 !~ ">" {c+=length($0);} END { print c; }\' > '
        + output_dir
        + bs
        + sample_name
        + "_len.txt && sed -i '/^$/d' "
        + output_dir
        + bs
        + sample_name
        + "_len.txt"
    )
    cmds = [cmd1]
    for cmd in cmds:
        subprocess.call(cmd, shell=True)
    scaftig_length_file = output_dir + "/" + sample_name + "_len.txt"

    fasta_file = old_pathname
    newfasta_file = output_dir + "/" + sample_name + "_scaf.fa"
    scaf_name = sample_name

    newnames = open(scaftig_length_file, "r")
    scaftig_file = open(newfasta_file, "w")

    fasta = open(fasta_file, "r")
    for line in fasta:
        if line.startswith(">"):
            newname = newnames.readline().split("\t")
            scaftig_file.write(">" + scaf_name + "_revised_" + newname[0] + "  " + newname[1])
        else:
            scaftig_file.write(line)
    scaftig_file.close()
    pathname = os.path.abspath(newfasta_file)
    ##    print 'Pathname:   ', pathname

    return scaftig_file


def header_replacer_F():
    global scaftig_file, pathname
    print("Manipulating assembly")
    ##    fasta_file = scaftig_file_unprocessed
    ##    newnames_file = scaftig_length_file
    newfasta_file = output_dir + "/" + sample_name + "_scaf.fa"
    scaf_name = sample_name
    fasta = open(scaftig_file_unprocessed, "r")
    newnames = open(scaftig_length_file, "r")
    scaftig_file = open(newfasta_file, "w")

    for line in fasta:
        if line.startswith(">"):
            newname = newnames.readline().split("\t")
            scaftig_file.write(">" + scaf_name + "_revised_" + newname[0] + "  " + newname[1])
        else:
            scaftig_file.write(line)
    pathname = os.path.abspath(newfasta_file)
    ##    print 'Pathname:   ', pathname
    ##    global pathname

    return scaftig_file


def gene_prediction(scaftig_file, sample_name):
    global gene_file
    print("gene prediction")
    scaftig_file = output_dir + "/" + sample_name + "_scaf.fa"
    cmd7 = (
        "prodigal -q -i "
        + scaftig_file
        + " -d "
        + output_dir
        + "/"
        + sample_name
        + "_seq -f gff -p meta -s "
        + output_dir
        + "/"
        + sample_name
        + "_potential_genes -o "
        + output_dir
        + "/"
        + sample_name
        + " > "
        + output_dir
        + "/gene_prediction.log"
    )
    cmds = [cmd7]
    for cmd in cmds:
        subprocess.call(cmd, shell=True)
    Prodigal_file = output_dir + "/" + sample_name
    Scaf_name, Gene_st, Gene_en, Gene_len, Strand, ID_info, Scaf_gene = [], [], [], [], [], [], []
    fasta = open(Prodigal_file, "r")
    gene_file = output_dir + "/" + sample_name + "_MC"

    ##this will convert the .gff prodigal output file to .tab file with information of our interest
    for line in fasta:
        if line.startswith("#"):
            continue
        else:
            temp_line = line.split("\t")
            G_inf = temp_line[8].split(";")
            G_ID = G_inf[0].split("_")
            scf_G = temp_line[0] + "_gene" + G_ID[1]
            Scaf_name.append(temp_line[0])
            Gene_st.append(temp_line[3])
            Gene_en.append(temp_line[4])
            Strand.append(temp_line[6])
            ID_info.append(G_ID[1])
            Scaf_gene.append(scf_G)

    Gene_st = [int(item) if item.isdigit() else item for item in Gene_st]
    Gene_en = [int(item) if item.isdigit() else item for item in Gene_en]

    for i in range(len(Gene_st)):
        length = Gene_en[i] - Gene_st[i]
        Gene_len.append(length)

    ##    Check if all the columns are of same length. If not then it means that one of the entries are missing.
    if (
        len(Scaf_name)
        == len(Scaf_gene) & len(Scaf_gene)
        == len(Gene_st) & len(Gene_st)
        == len(Gene_en) & len(Gene_en)
        == len(Gene_len) & len(Gene_len)
        == len(Strand)
    ):
        pass
    else:
        sys.exit(1)
        print("Number of entries for each gene are not equal, please correct the data and run again")

    ##    saving all lists to file
    zoo = zip(Scaf_gene, Strand, Gene_st, Gene_en, Gene_len)
    with open(gene_file, "wb") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerows(zoo)
    out_file.close()
    gene_pred_tool = ""
    gene_pred_tool == "2"
    print("gene prediction completed")
    return (gene_file, scaftig_file)


def gff2tab(sample_name):
    global gene_file

    Prodigal_file = gene_file
    Scaf_name, Gene_st, Gene_en, Gene_len, Strand, ID_info, Scaf_gene = [], [], [], [], [], [], []

    fasta = open(Prodigal_file, "r")
    gene_file = output_dir + "/" + sample_name + "_MC"

    ##this will convert the .gff prodigal output file to .tab file with information of our interest
    for line in fasta:
        if line.startswith("#"):
            continue
        else:
            temp_line = line.split("\t")
            G_inf = temp_line[8].split(";")
            G_ID = G_inf[0].split("_")
            scf_G = temp_line[0] + "_gene" + G_ID[1]
            Scaf_name.append(temp_line[0])
            Gene_st.append(temp_line[3])
            Gene_en.append(temp_line[4])
            Strand.append(temp_line[6])
            ID_info.append(G_ID[1])
            Scaf_gene.append(sample_name + "_revised_" + scf_G)

    Gene_st = [int(item) if item.isdigit() else item for item in Gene_st]
    Gene_en = [int(item) if item.isdigit() else item for item in Gene_en]

    for i in range(len(Gene_st)):
        length = Gene_en[i] - Gene_st[i]
        Gene_len.append(length)

    ##    Check if all the columns are of same length. If not then it means that one of the entries are missing.
    if (
        len(Scaf_name)
        == len(Scaf_gene) & len(Scaf_gene)
        == len(Gene_st) & len(Gene_st)
        == len(Gene_en) & len(Gene_en)
        == len(Gene_len) & len(Gene_len)
        == len(Strand)
    ):
        pass
    else:
        sys.exit(1)
        print("Number of entries for each gene are not equal, please correct the data and run again")

    ##    saving all lists to file
    zoo = zip(Scaf_gene, Strand, Gene_st, Gene_en, Gene_len)
    with open(gene_file, "w") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerows(zoo)
    out_file.close()
    gene_pred_tool = ""
    gene_pred_tool == "2"
    ##    pathname = scaftig_file

    return (gene_file, pathname)


def data_extraction(gene_file, gene_pred_tool):
    print("Gene data extraction")
    gene_name_tmp = []
    gene_name = []
    scaftig_no = []
    scaftig_name_tmp = []
    scaftig_name = []
    file_name = []
    scaftig_number = []
    gene_no = []
    just_genename = []
    scaf_name4dict = []
    scaf_name4dict_tmp = []
    ## extracting columns from the file
    gene_name_scaf = [x[0] for x in csv.reader(open(gene_file, "r"), delimiter="\t")]
    strand = [x[1] for x in csv.reader(open(gene_file, "r"), delimiter="\t")]
    gene_st = [x[2] for x in csv.reader(open(gene_file, "r"), delimiter="\t")]
    gene_end = [x[3] for x in csv.reader(open(gene_file, "r"), delimiter="\t")]

    #####################   metagenemark   #####################
    if gene_pred_tool == "1":
        print("METAGENEMARK")
        for i1 in gene_name_scaf:
            tmp = re.compile(
                "(?P<file_name>.*?)_*(?P<rev>revised)_*(?P<scaftig_name>[a-zA-Z0-9]*)_*(?P<scaftig_no>[0-9_]*)_(?P<gene_name>gene)(?P<number>\d+)"
            )
            match1 = tmp.search(i1)
            if match1:
                file_name.append(match1.group("file_name")), scaftig_name_tmp.append(
                    match1.group("file_name", "rev", "scaftig_name", "scaftig_no")
                ), scaf_name4dict_tmp.append(match1.group("file_name", "rev", "scaftig_name", "scaftig_no"))
                scaftig_no.append(match1.group("scaftig_no")), gene_name_tmp.append(
                    match1.group("gene_name", "scaftig_no", "number")
                ), just_genename.append(match1.group("gene_name")), gene_no.append(match1.group("number"))
        for i4 in scaftig_name_tmp:
            scaftig_name.append("_".join(i4))
        for i5 in gene_name_tmp:
            gene_name.append("_".join(i5))
        for i6 in scaf_name4dict_tmp:
            scaf_name4dict.append("_".join(i6))

    ##                                                  #####################   Prodigal    #####################
    elif gene_pred_tool == "2":
        print("PRODIGAL")
        for i2 in gene_name_scaf:
            tmp1 = re.compile(
                ##"(?P<file_name>.*?)_*(?P<rev>revised)_*(?P<scaftig_name>[a-zA-Z0-9]*)_(?P<gene_name>gene)(?P<scaftig_no>[0-9_]*)_(?P<gene_no>[0-9_]*)"
                "(?P<file_name>.*?)_*(?P<rev>revised)_*(?P<scaftig_name>[a-zA-Z0-9]*)_*(?P<scaftig_no>[0-9_]*)_(?P<gene_name>gene)(?P<gene_no>\d+)"
            )
            match2 = tmp1.search(i2)
            if match2:
                file_name.append(match2.group("file_name")), scaftig_name_tmp.append(
                    match2.group("file_name", "rev", "scaftig_name", "scaftig_no")
                ), scaf_name4dict_tmp.append(match2.group("file_name", "rev", "scaftig_name", "scaftig_no"))
                scaftig_no.append(match2.group("scaftig_no")), gene_name_tmp.append(
                    match2.group("gene_name", "scaftig_no", "gene_no")
                ), just_genename.append(match2.group("gene_name")), gene_no.append(match2.group("gene_no"))
        ##        print 'scaftig_name_tmp', scaftig_name_tmp[:10]
        for i6 in scaftig_name_tmp:
            scaftig_name.append("_".join(i6))
        for i7 in gene_name_tmp:
            gene_name.append("_".join(i7))
        for i8 in scaf_name4dict_tmp:
            scaf_name4dict.append("_".join(i8))

            ##converting str to int
    gene_st = map(int, gene_st)
    gene_end = map(int, gene_end)
    ##removing empty entries
    gene_name = [none1 for none1 in gene_name if none1 != ""]
    strand = [none2 for none2 in strand if none2 != ""]
    gene_st = [none3 for none3 in gene_st if none3 != ""]
    gene_end = [none4 for none4 in gene_end if none4 != ""]

    gene_name_scaf = [none1 for none1 in gene_name_scaf if none1 != ""]
    strand = [none2 for none2 in strand if none2 != ""]
    file_name = [none3 for none3 in file_name if none3 != ""]
    scaftig_name = [none4 for none4 in scaftig_name if none4 != ""]
    scaf_name4dict = [none4 for none4 in scaf_name4dict if none4 != ""]

    return gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict


def seq_info(gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict):
    print("Extracting sequence information...")
    ##    global pathname
    with open(pathname, mode="r", buffering=-1) as scaftig_file1:
        scaf_line = []
        scaf_len_only = []
        for SF in scaftig_file1:
            if SF.startswith(">"):
                scaf_line.append(SF[1:-1])
        scaffold_name_len = []
        for SL in scaf_line:
            scaffold_name_len.append(SL.split("length="))
        scaffold_name_length = [[SNL.strip() for SNL in inner] for inner in scaffold_name_len]
        for SLO in range(len(scaffold_name_length)):
            scaf_len_only.append(scaffold_name_length[SLO][1])
        scaf_len_only = [int(item) if item.isdigit() else item for item in scaf_len_only]
        scaf_name_len_dict = dict(scaffold_name_length)

        for key in scaf_name_len_dict:
            scaf_name_len_dict[key] = int(scaf_name_len_dict[key])
        ##        print scaf_name_len_dict, len(scaf_name_len_dict)

        if (
            len(scaftig_name)
            == len(gene_name) & len(gene_name)
            == len(strand) & len(strand)
            == len(gene_st) & len(gene_st)
            == len(gene_end) & len(gene_end)
            == len(scaf_name4dict)
        ):
            pass
            print("seq_info___Saving...")
            zoo = zip(scaftig_name, gene_name, strand, gene_st, gene_end)
            with open(output_dir + "/Gene_file_edited.tab", "w+") as out_file:
                writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
                writer.writerow(["scaftig_name", "gene_name", "strand", "gene_st", "gene_end"])
                writer.writerows(zoo)
            out_file.close()
        else:
            print("Number of entries for each gene are not equal, please correct the data and run again")
            exit()
    print("seqeuence info extraction___completed")
    return gene_name, strand, gene_st, gene_end, scaf_name_len_dict


def seq_info2(gene_name_scaf, strand, gene_st, gene_end, file_name, scaftig_name, gene_name, scaf_name4dict):
    print("Extracting sequence information...")
    ##    global pathname
    with open(pathname, mode="r", buffering=-1) as scaftig_file1:
        scaf_line = []
        scaf_len_only = []
        for SF in scaftig_file1:
            if SF.startswith(">"):
                scaf_line.append(SF[1:-1])
        scaffold_name_len = []
        for SL in scaf_line:
            scaffold_name_len.append(SL.split("length="))
        scaffold_name_length = [[SNL.strip() for SNL in inner] for inner in scaffold_name_len]
        for SLO in range(len(scaffold_name_length)):
            scaf_len_only.append(scaffold_name_length[SLO][1])
        scaf_len_only = [int(item) if item.isdigit() else item for item in scaf_len_only]
        scaf_name_len_dict = dict(scaffold_name_length)
        for key in scaf_name_len_dict:
            scaf_name_len_dict[key] = int(scaf_name_len_dict[key])

        if (
            len(scaftig_name)
            == len(gene_name) & len(gene_name)
            == len(strand) & len(strand)
            == len(gene_st) & len(gene_st)
            == len(gene_end) & len(gene_end)
            == len(scaf_name4dict)
        ):
            pass
            print("Seq_info___Saving...")
            zoo = zip(scaftig_name, gene_name, strand, gene_st, gene_end)
            with open(output_dir + "/Gene_file_edited.tab", "w+") as out_file:
                writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
                writer.writerow(["scaftig_name", "gene_name", "strand", "gene_st", "gene_end"])
                writer.writerows(zoo)
            out_file.close()
        else:
            print("Number of entries for each gene are not equal, please correct the data and run again")
            exit()
    print("seqeuence info extraction___completed")
    return gene_name, strand, gene_st, gene_end, scaf_name_len_dict


def upstream_coordinates_extraction(scaftig_name, gene_name, strand, gene_st, gene_end, scaf_name4dict):
    print("Calculating available upstream region...")
    upstream_seq_len = []
    scaftig_name_new = []
    gene_st_new = []
    gene_end_new = []
    strand_new = []
    gene_name_new = []
    new_scaf_name = []
    scaf_name4dict_new = []

    for i in range(len(scaftig_name)):
        if i < len(scaftig_name) - 1:
            if scaftig_name[i] != scaftig_name[i - 1] and scaftig_name[i] == scaftig_name[i + 1]:
                ups1 = gene_st[i] - 1
                upstream_seq_len.append(ups1)
                scaftig_name_new.append(scaftig_name[i])
                gene_name_new.append(gene_name[i])
                strand_new.append(strand[i])
                gene_st_new.append(gene_st[i])
                gene_end_new.append(gene_end[i])
                scaf_name4dict_new.append(scaf_name4dict[i])
            elif scaftig_name[i] == scaftig_name[i - 1]:
                ups1 = gene_st[i] - gene_end[i - 1]
                upstream_seq_len.append(ups1)
                scaftig_name_new.append(scaftig_name[i])
                gene_name_new.append(gene_name[i])
                strand_new.append(strand[i])
                gene_st_new.append(gene_st[i])
                gene_end_new.append(gene_end[i])
                scaf_name4dict_new.append(scaf_name4dict[i])
            elif scaftig_name[i] != scaftig_name[i + 1] and scaftig_name[i] == scaftig_name[i - 1]:
                ups1 = gene_st[i] - gene_end[i - 1]
                upstream_seq_len.append(ups1)
                scaftig_name_new.append(scaftig_name[i])
                gene_name_new.append(gene_name[i])
                strand_new.append(strand[i])
                gene_st_new.append(gene_st[i])
                gene_end_new.append(gene_end[i])
                scaf_name4dict_new.append(scaf_name4dict[i])
            elif scaftig_name[i] != scaftig_name[i + 1] and scaftig_name[i] != scaftig_name[i - 1]:
                continue
        if i == len(scaftig_name) - 1 and scaftig_name[i] == scaftig_name[i - 1]:
            ups1 = gene_st[i] - gene_end[i - 1]
            upstream_seq_len.append(ups1)
            scaftig_name_new.append(scaftig_name[i])
            gene_name_new.append(gene_name[i])
            strand_new.append(strand[i])
            gene_st_new.append(gene_st[i])
            gene_end_new.append(gene_end[i])
            scaf_name4dict_new.append(scaf_name4dict[i])
    upstream_seq_len = list(map(int, upstream_seq_len))
    gene_st_new = list(map(int, gene_st_new))

    print("Extracting upstream coordinates")
    upstream_st_pos = []
    upstream_en_pos = []
    for i in range(len(scaftig_name_new)):
        if i < len(scaftig_name_new) - 1:
            if i == 0:
                if upstream_seq_len[i] < 2:
                    up_st = 0
                    up_en = 0
                    upstream_st_pos.append(up_st)
                    upstream_en_pos.append(up_en)
                else:
                    up_st = gene_st_new[i] - upstream_seq_len[i]
                    upstream_st_pos.append(up_st)
                    up_en = gene_st_new[i] - 1
                    upstream_en_pos.append(up_en)
            elif i > 0 and i != len(scaftig_name_new) - 1:
                if upstream_seq_len[i] < 2:
                    up_st = 0
                    up_en = 0
                    upstream_st_pos.append(up_st)
                    upstream_en_pos.append(up_en)
                else:
                    if scaftig_name_new[i] != scaftig_name_new[i - 1]:
                        up_st = gene_st_new[i] - upstream_seq_len[i]
                        upstream_st_pos.append(up_st)
                        up_en = gene_st_new[i] - 1
                        upstream_en_pos.append(up_en)
                    elif scaftig_name_new[i] == scaftig_name_new[i + 1]:
                        up_st = gene_end_new[i - 1] + 1
                        upstream_st_pos.append(up_st)
                        up_en = gene_st_new[i] - 1
                        upstream_en_pos.append(up_en)
                    elif scaftig_name_new[i] != scaftig_name_new[i + 1]:
                        up_st = gene_end_new[i - 1] + 1
                        upstream_st_pos.append(up_st)
                        up_en = gene_st_new[i] - 1
                        upstream_en_pos.append(up_en)
        else:
            if scaftig_name_new[i] == scaftig_name_new[i - 1]:
                up_st = gene_end_new[i - 1] + 1
                upstream_st_pos.append(up_st)
                up_en = gene_st_new[i] - 1
                upstream_en_pos.append(up_en)
            else:
                up_st = 0
                up_en = 0
                upstream_st_pos.append(up_st)
                upstream_en_pos.append(up_en)
    return (
        upstream_seq_len,
        scaftig_name_new,
        gene_name_new,
        strand_new,
        gene_st_new,
        gene_end_new,
        scaf_name4dict_new,
        upstream_st_pos,
        upstream_en_pos,
    )


def scaftig_length(scaf_name4dict, scaf_name_len_dict):
    print("Extracting scaftig length information")
    leng_scaftig = []
    for i in range(len(scaf_name4dict)):
        scaf_len1 = scaf_name_len_dict.get(scaf_name4dict[i])
        leng_scaftig.append(scaf_len1)
    return leng_scaftig


def downstream_coordinate_region(scaftig_name, gene_st, gene_end, leng_scaftig, scaftig_name_new, gene_end_new):
    print("Calculating available downstream region...")
    downstream_seq_len = []
    downstream_st_pos = []
    downstream_en_pos = []

    for i in range(len(scaftig_name)):
        if i < len(scaftig_name) - 1:
            if scaftig_name[i] == scaftig_name[i + 1]:
                dss = gene_st[i + 1] - gene_end[i]
                downstream_seq_len.append(dss)
            elif i != 0 and scaftig_name[i] != scaftig_name[i + 1] and scaftig_name[i] == scaftig_name[i - 1]:
                dss = leng_scaftig[i] - gene_end[i]
                downstream_seq_len.append(dss)
            elif i == len(scaftig_name) and scaftig_name[i] == scaftig_name[i - 1]:
                dss = leng_scaftig[i] - gene_end[i]
                downstream_seq_len.append(dss)
                print(leng_scaftig[i] - gene_end[i])
            elif scaftig_name[i] != scaftig_name[i + 1] and scaftig_name[i] != scaftig_name[i - 1]:
                continue

        if i == len(scaftig_name) - 1 and scaftig_name[i] == scaftig_name[i - 1]:
            dss = leng_scaftig[i] - gene_end[i]
            downstream_seq_len.append(dss)
    print("Calculating downstream coordinates...")

    for i in range(len(scaftig_name_new)):
        if i < len(scaftig_name_new) - 1:
            if i == 0:
                if downstream_seq_len[i] < 2:
                    ds_st = 0
                    ds_en = 0
                    downstream_st_pos.append(ds_st)
                    downstream_en_pos.append(ds_en)
                else:
                    ds_st = gene_end_new[i] + 1
                    ds_en = (gene_end_new[i] + downstream_seq_len[i]) - 1
                    downstream_st_pos.append(ds_st)
                    downstream_en_pos.append(ds_en)
            elif i > 0 and i != len(scaftig_name_new) - 1:
                if downstream_seq_len[i] < 2:
                    ds_st = 0
                    ds_en = 0
                    downstream_st_pos.append(ds_st)
                    downstream_en_pos.append(ds_en)
                else:
                    if scaftig_name_new[i] != scaftig_name_new[i - 1]:
                        ds_st = gene_end_new[i] + 1
                        ds_en = (gene_end_new[i] + downstream_seq_len[i]) - 1
                        downstream_st_pos.append(ds_st)
                        downstream_en_pos.append(ds_en)
                    elif scaftig_name_new[i] == scaftig_name_new[i + 1]:
                        ds_st = gene_end_new[i] + 1
                        ds_en = (gene_end_new[i] + downstream_seq_len[i]) - 1
                        downstream_st_pos.append(ds_st)
                        downstream_en_pos.append(ds_en)
                    elif scaftig_name_new[i] != scaftig_name_new[i + 1]:
                        ds_st = gene_end_new[i] + 1
                        ds_en = (gene_end_new[i] + downstream_seq_len[i]) - 1
                        downstream_st_pos.append(ds_st)
                        downstream_en_pos.append(ds_en)
        else:
            if scaftig_name_new[i] == scaftig_name_new[i - 1]:
                ds_st = gene_end_new[i] + 1
                ds_en = gene_end_new[i] + downstream_seq_len[i]
                downstream_st_pos.append(ds_st)
                downstream_en_pos.append(ds_en)
            else:
                ds_st = 0
                ds_en = 0
                downstream_st_pos.append(ds_st)
                downstream_en_pos.append(ds_en)
    return downstream_seq_len, downstream_st_pos, downstream_en_pos


def saving_data(
    scaftig_name_new, gene_name_new, strand_new, gene_st_new, gene_end_new, upstream_seq_len, downstream_seq_len
):
    print("Upstream & Downstream Length Information___Saving")
    if (
        len(scaftig_name_new)
        == len(gene_name_new) & len(gene_name_new)
        == len(strand_new) & len(strand_new)
        == len(gene_st_new) & len(gene_st_new)
        == len(gene_end_new) & len(gene_end_new)
        == len(upstream_seq_len) & len(upstream_seq_len)
        == len(downstream_seq_len)
    ):
        pass
        zoo = zip(
            scaftig_name_new, gene_name_new, strand_new, gene_st_new, gene_end_new, upstream_seq_len, downstream_seq_len
        )
        with open(output_dir + "/Gene_file_UDlen.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerow(["scaftig_name", "gene_name", "strand", "gene_st", "gene_end", "upst_len", "dnst_len"])
            writer.writerows(zoo)
        out_file.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()
    return


def saving_data2(
    scaftig_name_new,
    gene_name_new,
    strand_new,
    gene_st_new,
    gene_end_new,
    upstream_seq_len,
    downstream_seq_len,
    upstream_st_pos,
    upstream_en_pos,
    downstream_st_pos,
    downstream_en_pos,
):
    print("All Genes Info___Saving")
    ##SAVING ALL DATA TO A NEWFILE

    if (
        len(scaftig_name_new) == len(gene_name_new)
        and len(gene_name_new) == len(strand_new)
        and len(strand_new) == len(gene_st_new)
        and len(gene_st_new) == len(gene_end_new)
        and len(gene_end_new) == len(upstream_seq_len)
        and len(upstream_seq_len) == len(downstream_seq_len)
        and len(downstream_seq_len) == len(upstream_st_pos)
        and len(upstream_st_pos) == len(upstream_en_pos)
        and len(upstream_en_pos) == len(downstream_st_pos)
        and len(downstream_st_pos) == len(downstream_en_pos)
    ):
        zoo = zip(
            scaftig_name_new,
            gene_name_new,
            strand_new,
            gene_st_new,
            gene_end_new,
            upstream_seq_len,
            upstream_st_pos,
            upstream_en_pos,
            downstream_seq_len,
            downstream_st_pos,
            downstream_en_pos,
        )
        with open(output_dir + "/Gene_file_detailed.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Scaf_name",
                    "Gene_name",
                    "Strand",
                    "Gene_st",
                    "Gene_en",
                    "UpSt_Reg",
                    "UP_St",
                    "uP_En",
                    "DnSt_reg",
                    "Dn_St",
                    "Dn_En",
                ]
            )
            writer.writerows(zoo)
        out_file.close()
    else:
        print("data not correct")
        exit
    return


def UPS_DSS_Slicing(
    downstream_seq_len,
    downstream_st_pos,
    downstream_en_pos,
    scaftig_name_new,
    gene_name_new,
    upstream_seq_len,
    upstream_st_pos,
    upstream_en_pos,
):
    ###### MANIPULATING LENGHT OF COORDINATES
    print("UPS_DSS_Slicing")
    newUPS_st_pos = []
    newUPS_en_pos = []
    newDSS_st_pos = []
    newDSS_en_pos = []
    gene_name4cord_DSS = []
    scaf_name4cord_DSS = []
    scaf_name4cord_UPS = []
    gene_name4cord_UPS = []

    for i in range(len(scaftig_name_new)):
        if downstream_seq_len[i] >= 15:
            if downstream_seq_len[i] <= 700:
                dss_st_new = downstream_st_pos[i]
                dss_en_new = downstream_en_pos[i]
                newDSS_st_pos.append(dss_st_new)
                newDSS_en_pos.append(dss_en_new)
                scaf_name4cord_DSS.append(scaftig_name_new[i])
                gene_name4cord_DSS.append(gene_name_new[i])
            elif downstream_seq_len[i] > 700:
                dss_st_new = downstream_st_pos[i]
                dss_en_new = downstream_st_pos[i] + 700
                newDSS_st_pos.append(dss_st_new)
                newDSS_en_pos.append(dss_en_new)
                scaf_name4cord_DSS.append(scaftig_name_new[i])
                gene_name4cord_DSS.append(gene_name_new[i])
        elif downstream_seq_len[i] < 15:
            dss_st_new = "short_dss"
            dss_en_new = "short_dss"
            newDSS_st_pos.append(dss_st_new)
            newDSS_en_pos.append(dss_en_new)
            scaf_name4cord_DSS.append(scaftig_name_new[i])
            gene_name4cord_DSS.append(gene_name_new[i])
        if upstream_seq_len[i] >= 15:
            if upstream_seq_len[i] <= 700:
                ups_st_new = upstream_st_pos[i]
                ups_en_new = upstream_en_pos[i]
                newUPS_st_pos.append(ups_st_new)
                newUPS_en_pos.append(ups_en_new)
                scaf_name4cord_UPS.append(scaftig_name_new[i])
                gene_name4cord_UPS.append(gene_name_new[i])
            elif upstream_seq_len[i] > 700:
                ups_st_new = upstream_en_pos[i] - 700
                ups_en_new = upstream_en_pos[i]
                newUPS_st_pos.append(ups_st_new)
                newUPS_en_pos.append(ups_en_new)
                scaf_name4cord_UPS.append(scaftig_name_new[i])
                gene_name4cord_UPS.append(gene_name_new[i])
        elif upstream_seq_len[i] < 15:
            ups_st_new = "short_ups"
            ups_en_new = "short_ups"
            newUPS_st_pos.append(ups_st_new)
            newUPS_en_pos.append(ups_en_new)
            scaf_name4cord_UPS.append(scaftig_name_new[i])
            gene_name4cord_UPS.append(gene_name_new[i])

    return (
        newDSS_st_pos,
        newDSS_en_pos,
        newUPS_st_pos,
        newUPS_en_pos,
        gene_name4cord_DSS,
        scaf_name4cord_DSS,
        scaf_name4cord_UPS,
        gene_name4cord_UPS,
    )


def Data_saving3(
    newDSS_st_pos,
    newDSS_en_pos,
    newUPS_st_pos,
    newUPS_en_pos,
    gene_name4cord_DSS,
    scaf_name4cord_DSS,
    scaf_name4cord_UPS,
    gene_name4cord_UPS,
    scaftig_name_new,
    gene_name_new,
    strand_new,
    gene_st_new,
    gene_end_new,
):
    print("All Gene Coordinates___Saving")

    ####SAVING COORDINATES TO FILES
    if (
        len(newUPS_st_pos) == len(newUPS_en_pos)
        and len(newUPS_en_pos) == len(newDSS_st_pos)
        and len(newDSS_st_pos) == len(strand_new)
        and len(strand_new) == len(newDSS_en_pos)
        and len(newDSS_en_pos) == len(scaf_name4cord_UPS)
        and len(scaf_name4cord_UPS) == len(gene_name4cord_UPS)
        and len(gene_name4cord_UPS) == len(gene_name4cord_DSS)
        and len(gene_name4cord_DSS) == len(scaf_name4cord_DSS)
    ):
        zoo = zip(
            scaftig_name_new,
            gene_name_new,
            strand_new,
            gene_st_new,
            gene_end_new,
            newUPS_st_pos,
            newUPS_en_pos,
            newDSS_st_pos,
            newDSS_en_pos,
        )
        with open(output_dir + "/Gene_coord.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["Scaf_name", "Gene_name", "Strand", "Gene_st", "Gene_en", "UP_St", "uP_En", "Dn_St", "Dn_En"]
            )
            writer.writerows(zoo)
        out_file.close()
    else:
        print("data not correct")
        exit
    #### saving UPS and DSS in seperate files
    if (
        len(scaftig_name_new) == len(gene_name_new)
        and len(gene_name_new) == len(newUPS_st_pos)
        and len(newUPS_st_pos) == len(newUPS_en_pos)
        and len(newUPS_en_pos)
    ):
        zoo = zip(scaftig_name_new, gene_name_new, newUPS_st_pos, newUPS_en_pos)
        with open(output_dir + "/Gene_ups_Coordinates.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerow(["Scaf_name", "Gene_name", "UP_St", "uP_En"])
            writer.writerows(zoo)
        out_file.close()
        upstream_coordinate_file = output_dir + "/Gene_ups_Coordinates.tab"
    if (
        len(scaftig_name_new) == len(gene_name_new)
        and len(gene_name_new) == len(newDSS_st_pos)
        and len(newDSS_st_pos) == len(newDSS_en_pos)
        and len(newDSS_en_pos)
    ):
        zoo = zip(scaftig_name_new, gene_name_new, newDSS_st_pos, newDSS_en_pos)
        with open(output_dir + "/Gene_dss_Coordinates.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerow(["Scaf_name", "Gene_name", "DSS_St", "DSS_En"])
            writer.writerows(zoo)
        out_file.close()
        downstream_coordinate_file = output_dir + "/Gene_dss_Coordinates.tab"
    return downstream_coordinate_file, upstream_coordinate_file


def getsouce():
    print("Sequence Extraction")
    sourcedic = {}
    sourcefile = pathname
    with open(sourcefile) as f:
        sourcestring = f.read().strip(">")
    removelist = sourcestring.split(">")
    for removestri in removelist:
        genlist = removestri.strip().split("\n")
        genNamestr = genlist[0].split("  ")
        genlist.remove(genlist[0])
        sourcedic[genNamestr[0]] = "".join(genlist)
    return sourcedic


def getsouce_op(pathname):
    print("Sequence Extraction___Processing")
    sourcedic = {}
    sourcefile = pathname
    ##    sourcefile = open(scaf_file, 'r')
    ##	sourcefile = scaftig_file1
    with open(sourcefile) as f:
        sourcestring = f.read().strip(">")
    removelist = sourcestring.split(">")
    for removestri in removelist:
        genlist = removestri.strip().split("\n")
        genNamestr = genlist[0].split("  ")
        genlist.remove(genlist[0])
        sourcedic[genNamestr[0]] = "".join(genlist)
    ##    print 'getsouce completed'
    return sourcedic


def getgenstring_ups(gendic, genefile):
    print("Upstream___Sequence Extraction")
    ##    print 'Upstream sequence extraction'
    with open(output_dir + "/Upstream_seq.fasta", "w+") as out_ups:
        with open(genefile) as f:
            geneslist = f.read().strip().split("\n")
        geneslist.remove(geneslist[0])
        for index in geneslist:
            outstr = index.split("\t")
            namelist = index.split("\t")
            if namelist[2] != "short_ups":
                name_ups = ">" + namelist[0] + "," + namelist[1] + " " + namelist[2] + "\t" + namelist[-1] + "\n"
                ups_str = gendic[outstr[0]][(int(outstr[2]) - 1) : (int(outstr[3]) - 1)] + "\n"
                out_ups.write(name_ups)
                out_ups.write(ups_str)
    upstream_seq_file = open(output_dir + "/Upstream_seq.fasta", "r")
    ##    print 'getgenstring_ups completed'
    return upstream_seq_file


def getgenstring_dss(gendic, genefile):
    print("Downstream___Sequence Extraction")
    with open(output_dir + "/Downstream_seq.fasta", "w+") as out_dss:
        with open(genefile) as f:
            geneslist = f.read().strip().split("\n")
        geneslist.remove(geneslist[0])
        for index in geneslist:
            outstr = index.split("\t")
            namelist = index.split("\t")
            if namelist[2] != "short_dss":
                name_dss = ">" + namelist[0] + "," + namelist[1] + " " + namelist[2] + "\t" + namelist[-1] + "\n"
                dss_str = gendic[outstr[0]][(int(outstr[2]) - 1) : (int(outstr[3]) - 1)] + "\n"
                out_dss.write(name_dss)
                out_dss.write(dss_str)
    downstream_seq_file = open(output_dir + "/Downstream_seq.fasta", "r")
    ##    print 'Downstream sequence extraction start completed'
    return downstream_seq_file


def IGD_calc(scaftig_name_new, strand_new, gene_name_new, gene_st_new, gene_end_new):
    print("Calculating Intergenic Distance...")
    IGD, new_gene1_name, new_gene2_name, new_strand_1, new_strand_2, new_scaftig_no, new_scaftig_no, new_k, cls = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    (
        IGDclstr,
        new_gene1_clstr,
        new_gene2_clstr,
        new_strand_1clstr,
        new_strand_2clstr,
        new_scaftig_no_clstr,
        new_kclstr,
    ) = ([], [], [], [], [], [], [])
    k = 0
    c = 1
    ##print 'Gene name,\tStrand\tGene2_st\tGene1_end\tIGD\n'
    ##gene_st_new, gene_end_new,gene_name_new
    for i in range(len(scaftig_name_new[:-1])):
        if scaftig_name_new[i] == scaftig_name_new[i + 1]:
            if strand_new[i] == strand_new[i + 1]:
                a = gene_st_new[i + 1] - gene_end_new[i]
                IGD.append(a)
                new_gene1_name.append(gene_name_new[i])
                new_gene2_name.append(gene_name_new[i + 1])
                new_strand_1.append(strand_new[i])
                new_strand_2.append(strand_new[i + 1])
                new_scaftig_no.append(scaftig_name_new[i])
                k = i + 1
                new_k.append(k)
        if scaftig_name_new[i] == scaftig_name_new[i + 1]:
            if strand_new[i] == strand_new[i + 1]:
                if (
                    gene_st_new[i + 1] - gene_end_new[i] <= 600
                ):  ##                a = gene_st_new[i+1] - gene_end_new[i]
                    IGDclstr.append(a)
                    new_gene1_clstr.append(gene_name_new[i])
                    new_gene2_clstr.append(gene_name_new[i + 1])
                    new_strand_1clstr.append(strand_new[i])
                    new_strand_2clstr.append(strand_new[i + 1])
                    new_scaftig_no_clstr.append(scaftig_name_new[i])
                    k = i + 1
                    new_kclstr.append(k)

    ##            print k, scaftig_name_new[i], gene_name_new[i], ',',gene_name_new[i+1], '\t',strand_new[i],',',strand_new[i+1],'\t', gene_st_new[i+1],'\t', gene_end_new[i],'\t',a
    ##    print "IGD calculation done."
    return (
        new_k,
        new_scaftig_no,
        new_gene1_name,
        new_gene2_name,
        new_strand_1,
        new_strand_2,
        IGD,
        IGDclstr,
        new_gene1_clstr,
        new_gene2_clstr,
        new_strand_1clstr,
        new_strand_2clstr,
        new_scaftig_no_clstr,
        new_kclstr,
    )


def Data_saving4(
    new_k,
    new_scaftig_no,
    new_gene1_name,
    new_gene2_name,
    new_strand_1,
    new_strand_2,
    IGD,
    IGDclstr,
    new_gene1_clstr,
    new_gene2_clstr,
    new_strand_1clstr,
    new_strand_2clstr,
    new_scaftig_no_clstr,
    new_kclstr,
):
    print("Inergenic Distance Information/ Proximon___Saving")
    ### writing IGD& DIRECTIONAL CLUSTERS TO FILE
    zoo1 = zip(new_k, new_scaftig_no, new_gene1_name, new_gene2_name, new_strand_1, new_strand_2, IGD)
    with open(output_dir + "/DIR_Proximons_temp.tab", "w+") as IGD_OUTFILE:
        writer_IGD = csv.writer(IGD_OUTFILE, delimiter="\t", lineterminator="\n")
        writer_IGD.writerows(zoo1)
    zoo1 = zip(new_k, new_scaftig_no, new_gene1_name, new_gene2_name, new_strand_1, new_strand_2, IGD)
    with open(output_dir + "/DIR_Proximons_hdr.tab", "w+") as IGD_OUTFILE:
        writer_IGD = csv.writer(IGD_OUTFILE, delimiter="\t", lineterminator="\n")
        writer_IGD.writerow(["#", "SC_nm", "Gene1,Gene2", "Strand1,Strand2", "IGD"])
        writer_IGD.writerows(zoo1)
    zoo1 = zip(
        new_kclstr,
        new_scaftig_no_clstr,
        new_gene1_clstr,
        new_gene2_clstr,
        new_strand_1clstr,
        new_strand_2clstr,
        IGDclstr,
    )
    with open(output_dir + "/Proximons_hdr.tab", "w+") as IGD_OUTFILE:
        writer_IGD = csv.writer(IGD_OUTFILE, delimiter="\t", lineterminator="\n")
        writer_IGD.writerow(["#", "SC_nm", "Gene1,Gene2", "Strand1,Strand2", "IGD"])
        writer_IGD.writerows(zoo1)
    ##    print 'Data saving 4 done'
    IGD_OUTFILE.close()


def promoter_prediction(upstream_seq_file):
    global NNPP2_path
    print("Promoter Prediction")

    try:
        os.remove(output_dir + "/promoter_prediction.txt")
    except OSError:
        pass
    if os.path.exists(output_dir + "/promoter_prediction.txt"):
        os.remove(output_dir + "/promoter_prediction.txt")
        print("Removing previous Promoter file")
        promoter_file = open(output_dir + "/promoter_prediction.txt", "w")
    else:
        pass
    ##
    ##    print 'Promoter Prediction...'
    ######Extracting promoter sequence one by one
    promoter_file = open(output_dir + "/promoter_prediction.txt", "w")
    ##    loops = 0
    for line in upstream_seq_file:
        line = line.rstrip()
        if line[0] == ">":
            header = line[0:]
        ##            loops += 1
        else:
            sequence = line[0:]
            with open(output_dir + "/tmp_seq.txt", "w") as outfile:
                outfile.write(header)
                outfile.write("\n")
                outfile.write(sequence)
                outfile.write("\n")
        cmd1 = "cd " + NNPP2_path
        cmd2 = (
            NNPP2_path
            + "/bin/fa2TDNNpred-PRO-static.linux -t 0.25 "
            + output_dir
            + "/tmp_seq.txt >> %s/promoter_prediction.txt" % output_dir
        )
        ##    cmd4 = 'rm /home/alisu/Software/Promoter_prediction/NNPP/NNPP2.2/test/tmp_seq.txt'
        cmds = [cmd1, cmd2]
        for cmd in cmds:
            subprocess.call(cmd, shell=True)
    promoter_file.close()
    ##    promoter_file_path = output_dir+"/promoter_prediction.txt"
    ##    print 'Total Entries processed for promoter prediction: ', loops
    ##    print 'promoter_prediction'
    return


def Promoter_file_parse():
    promoter_file_path2 = open(output_dir + "/promoter_prediction.txt", "r")
    print("Parsing promoter file...")
    name = []
    answer = []
    loops = 0
    for i in promoter_file_path2:
        i = i.rstrip()
        if i.startswith("Prediction for"):
            Name1 = i[15:]
            name.append(Name1)
        elif i.startswith("No Hits with this threshold"):
            NoProm = "No"
            answer.append(NoProm)
        elif i.startswith("Hit 1:"):
            YesProm = "Yes"
            loops += 1
            answer.append(YesProm)
        else:
            pass
    name = name[1:]
    answer = answer[1:]
    print("Total Promoters Predicter:  ", loops)
    print("Promoter Parsed ___Saving")
    if len(name) == len(answer):
        zoo = zip(name, answer)
        with open(output_dir + "/Promoter_1.tab", "w+") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        outfile.close()
    scf_nm = []
    prom_res = []
    gene_nm = []

    counter = 0
    promoter_parsed_file = open(output_dir + "/Promoter_1.tab", "r")
    for i in promoter_parsed_file:
        i = i.rstrip()
        counter += 1
        tmp1 = re.split(",| \(|:\\t", i)
        scf_nm.append(tmp1[0])
        gene_nm.append(tmp1[1])
        prom_res.append(tmp1[-1])
    ##    print 'Parsing promoter file completed ...'

    return scf_nm, gene_nm, prom_res, counter


##--------------------------------------------------------------------------


def Prom_IGD_Clustering():
    print("Clustering___Promoter-Intergenic Distance")
    ##    print 'File Parsing...'
    igd_file = output_dir + "/DIR_Proximons_temp.tab"

    prom_scf_gene, prom_res = [], []
    counter = 0
    promoter_parsed_file = open(output_dir + "/Promoter_1.tab", "r")

    for i in promoter_parsed_file:
        i = i.rstrip()
        counter += 1
        tmp1 = re.split(" |\):\t", i)
        prom_scf_gene.append(tmp1[0])
        prom_res.append(tmp1[-1])

    if len(prom_scf_gene) == len(prom_res):
        print("Saving...")
        zoo = zip(prom_scf_gene, prom_res)
        with open(output_dir + "/Promoter_edited.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        out_file.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()

    igd_newk = [x[0] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_scf_nm = [x[1] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_gene1 = [x[2] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_gene2 = [x[3] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_st1 = [x[4] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_st2 = [
        x[5]
        for x in csv.reader(
            open(igd_file, "r"),
            delimiter="\t",
        )
    ]
    IGD = [x[6] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]

    igd_scf_gene1 = []
    igd_scf_gene2 = []
    for i in range(len(igd_scf_nm)):
        tmp1 = igd_scf_nm[i] + "," + igd_gene1[i]
        tmp2 = igd_scf_nm[i] + "," + igd_gene2[i]
        igd_scf_gene1.append(tmp1)
        igd_scf_gene2.append(tmp2)

    print("Prom_IGD_Clustering____Promter search...")
    (
        princ_sc_nm1,
        princ_sc_nm2,
        res_prom1,
        res_prom2,
        newk1_igd,
        newk2_igd,
        scnm1_igd,
        scnm2_igd,
        gene1igd,
        gene2igd,
        st1_igd,
        st2_igd,
        IGD1_igd,
        IGD2_igd,
    ) = ([], [], [], [], [], [], [], [], [], [], [], [], [], [])
    no_ans = "NA"
    counter = 0
    counter2 = 0
    counter3 = 0
    counter4 = 0
    counter5 = 0

    for i, j in zip(igd_scf_gene1, igd_scf_gene2):
        if i in prom_scf_gene:
            indx1_igd = igd_scf_gene1.index(i)
            indx_prom1 = prom_scf_gene.index(i)
            counter += 1
            newk1_igd.append(igd_newk[indx1_igd]), princ_sc_nm1.append(i), scnm1_igd.append(
                igd_scf_nm[indx1_igd]
            ), gene1igd.append(igd_gene1[indx1_igd]), st1_igd.append(igd_st1[indx1_igd]),
            res_prom1.append(prom_res[indx_prom1]), IGD1_igd.append(IGD[indx1_igd])
        ##        print counter, igd_newk[indx1_igd], i, '\t',  igd_scf_nm[indx1_igd], igd_gene1[indx1_igd], igd_st1[indx1_igd], prom_res[indx_prom1]
        elif i not in prom_scf_gene:
            indx1_igd = igd_scf_gene1.index(i)
            counter += 1
            newk1_igd.append(igd_newk[indx1_igd]), princ_sc_nm1.append(i), scnm1_igd.append(
                igd_scf_nm[indx1_igd]
            ), gene1igd.append(igd_gene1[indx1_igd]), st1_igd.append(igd_st1[indx1_igd]), res_prom1.append(
                no_ans
            ), IGD1_igd.append(
                IGD[indx1_igd]
            )
        ##        print counter, igd_newk[indx1_igd], i, '\t', igd_scf_nm[indx1_igd], igd_gene1[indx1_igd], igd_st1[indx1_igd], no_ans
        if j in prom_scf_gene:
            indx2_igd = igd_scf_gene2.index(j)
            indx_prom2 = prom_scf_gene.index(j)
            counter += 1
            newk2_igd.append(igd_newk[indx2_igd]), princ_sc_nm2.append(j), scnm2_igd.append(
                igd_scf_nm[indx2_igd]
            ), gene2igd.append(igd_gene2[indx2_igd]), st2_igd.append(igd_st2[indx2_igd]), res_prom2.append(
                prom_res[indx_prom2]
            ), IGD2_igd.append(
                IGD[indx2_igd]
            )
        ##        print counter, igd_newk[indx2_igd], j, '\t', igd_scf_nm[indx2_igd], igd_gene2[indx2_igd], igd_st2[indx2_igd], prom_res[indx_prom2]
        elif j not in prom_scf_gene:
            counter += 1
            indx2_igd = igd_scf_gene2.index(j)
            newk2_igd.append(igd_newk[indx2_igd]), princ_sc_nm2.append(j), scnm2_igd.append(
                igd_scf_nm[indx2_igd]
            ), gene2igd.append(igd_gene2[indx2_igd]), st2_igd.append(igd_st2[indx2_igd]), res_prom2.append(
                no_ans
            ), IGD2_igd.append(
                IGD[indx2_igd]
            )
    ##        print counter, igd_newk[indx2_igd], j, '\t', igd_scf_nm[indx2_igd], igd_gene2[indx2_igd], igd_st2[indx2_igd], no_ans
    ##    print '\n\n\nPROCESSING AFTER PROMOTER PREDICTION'
    ##    for i in range(20):
    ##        print newk1_igd[i],newk2_igd[i], princ_sc_nm1[i],princ_sc_nm2[i],  scnm1_igd[i],scnm2_igd[i],  gene1igd[i],gene2igd[i],  st1_igd[i], st2_igd[i], res_prom1[i], res_prom2[i], IGD1_igd[i], IGD2_igd[i]

    ############            CHECK

    if len(igd_scf_gene1) == len(newk2_igd) and len(newk2_igd) == len(newk1_igd):
        for i in range(len(newk2_igd)):
            if (
                newk1_igd[i] == newk2_igd[i]
                and scnm1_igd[i] == scnm2_igd[i]
                and gene1igd[i] != gene2igd[i]
                and st1_igd[i] == st2_igd[i]
                and IGD1_igd[i] == IGD2_igd[i]
            ):
                pass
            else:
                "Values not equal"
    else:
        print("length not equal")

    ##    print '\n\n\nPROCESSING 2'
    ##    for i in range(20):
    ##        print newk1_igd[i],newk2_igd[i], scnm1_igd[i], scnm2_igd[i], gene1igd[i], gene2igd[i], st1_igd[i], st2_igd[i],IGD1_igd[i],IGD2_igd[i]

    igd_clstr = newk1_igd
    igd_scaffold_name = scnm1_igd
    igd_geneA = gene1igd
    igd_geneB = gene2igd
    intergenic_dist = IGD1_igd
    igd_stA = st1_igd
    igd_stB = st2_igd
    promA_res = res_prom1
    promB_res = res_prom2

    if (
        len(igd_clstr) == len(igd_scaffold_name)
        and len(igd_scaffold_name) == len(igd_geneA)
        and len(igd_geneA) == len(igd_geneB)
        and len(igd_geneB) == len(intergenic_dist)
        and len(intergenic_dist) == len(igd_stA)
        and len(igd_stA) == len(igd_stB)
        and len(igd_stB) == len(promA_res)
        and len(promA_res) == len(promB_res)
    ):
        pass
        print("Prom_IGD_Clustering___Saving...")
        zoo = zip(
            igd_clstr, igd_scaffold_name, igd_geneA, igd_geneB, intergenic_dist, igd_stA, igd_stB, promA_res, promB_res
        )
        with open(output_dir + "/IGD_PROM_RES.tab", "w+") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        outfile.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()

    ##    print '\n\n\nPROCESSING 2 IMPO'
    ##    for i in range(20):
    ##        print igd_clstr[i],igd_scaffold_name[i], igd_geneA[i], igd_geneB[i], intergenic_dist[i], igd_stA[i], igd_stB[i], promA_res[i],promB_res[i]

    countera = 0
    (
        new_igd_clstr1,
        new_igd_scaffold_name1,
        new_igd_geneA1,
        new_igd_geneB1,
        new_intergenic_dist1,
        new_igd_stA1,
        new_igd_stB1,
        new_promA_res1,
        new_promB_res1,
    ) = ([], [], [], [], [], [], [], [], [])
    for i in range(len(igd_clstr)):
        if promA_res[i] == "Yes" and promB_res[i] == "No" or promB_res[i] == "NA":
            new_igd_clstr1.append(igd_clstr[i]), new_igd_scaffold_name1.append(
                igd_scaffold_name[i]
            ), new_igd_geneA1.append(igd_geneA[i]), new_igd_geneB1.append(igd_geneB[i]),
            new_intergenic_dist1.append(intergenic_dist[i]), new_igd_stA1.append(igd_stA[i]), new_igd_stB1.append(
                igd_stB[i]
            ), new_promA_res1.append(promA_res[i]), new_promB_res1.append(promB_res[i])
        elif promA_res[i] == "No" or promA_res[i] == "NA" and promB_res[i] == "No" or promB_res[i] == "NA":
            new_igd_clstr1.append(igd_clstr[i]), new_igd_scaffold_name1.append(
                igd_scaffold_name[i]
            ), new_igd_geneA1.append(igd_geneA[i]), new_igd_geneB1.append(igd_geneB[i]),
            new_intergenic_dist1.append(intergenic_dist[i]), new_igd_stA1.append(igd_stA[i]), new_igd_stB1.append(
                igd_stB[i]
            ), new_promA_res1.append(promA_res[i]), new_promB_res1.append(promB_res[i])
        elif promB_res[i] == "Yes" and promA_res[i] == "No" or promA_res[i] == "NA":
            countera += 1
            pass
        else:
            countera += 1
            pass
    ##    print 'igd_clstr LEN', len(igd_clstr)
    ##    print 'new_igd_clstr LEN', len(new_igd_clstr1)
    ##    print 'promA and B len', len(new_promA_res1) , len(new_promB_res1)
    ##    print '\n\n\nPROCESSING 3 IMPO'
    ##    for i in range(20):
    ##        print new_igd_clstr1[i],new_igd_scaffold_name1[i], new_igd_geneA1[i], new_igd_geneB1[i], new_intergenic_dist1[i], new_igd_stA1[i], new_igd_stB1[i], new_promA_res1[i],new_promB_res1[i]
    (
        new_igd_clstr,
        new_igd_scaffold_name,
        new_igd_geneA,
        new_igd_geneB,
        new_intergenic_dist,
        new_igd_stA,
        new_igd_stB,
        new_promA_res,
        new_promB_res,
    ) = ([], [], [], [], [], [], [], [], [])
    for i in range(len(new_promA_res1)):
        ##        if new_promA_res1[i] == 'Yes' and new_promB_res1[i] == 'No' or new_promB_res1[i] == 'NA':
        ##            new_igd_clstr.append(new_igd_clstr1[i]), new_igd_scaffold_name.append(new_igd_scaffold_name1[i]), new_igd_geneA.append(new_igd_geneA1[i]), new_igd_geneB.append(new_igd_geneB1[i]), new_intergenic_dist.append(new_intergenic_dist1[i]), new_igd_stA.append(new_igd_stA[i]), new_igd_stB.append(new_igd_stB[i]), new_promA_res.append(new_promA_res[i]), new_promB_res.append(new_promB_res[i])
        if new_promA_res1[i] == "No" or new_promA_res1[i] == "NA" and new_promB_res1[i] == "Yes":
            pass
        else:
            new_igd_clstr.append(new_igd_clstr1[i]), new_igd_scaffold_name.append(
                new_igd_scaffold_name1[i]
            ), new_igd_geneA.append(new_igd_geneA1[i]), new_igd_geneB.append(
                new_igd_geneB1[i]
            ), new_intergenic_dist.append(
                new_intergenic_dist1[i]
            ), new_igd_stA.append(
                new_igd_stA1[i]
            ), new_igd_stB.append(
                new_igd_stB1[i]
            ), new_promA_res.append(
                new_promA_res1[i]
            ), new_promB_res.append(
                new_promB_res1[i]
            )
    ##            del new_igd_clstr[i]
    ##            del new_igd_scaffold_name[i]
    ##            del new_igd_geneA[i]
    ##            del new_igd_geneB[i]
    ##            del new_intergenic_dist[i]
    ##            del new_igd_stA[i]
    ##            del new_igd_stB[i]
    ##            del new_promA_res[i]
    ##            del new_promB_res[i]

    ##    print 'igd_clstr LEN', len(igd_clstr)
    ##    print 'new_igd_clstr1 LEN', len(new_igd_clstr1)
    ##    print 'new_igd_clstr LEN', len(new_igd_clstr)
    ##    print '\n\n\nPROCESSING 4 IMPO'
    ##    for i in range(20):
    ##        print new_igd_clstr[i],new_igd_scaffold_name[i], new_igd_geneA[i], new_igd_geneB[i], new_intergenic_dist[i], new_igd_stA[i], new_igd_stB[i], new_promA_res[i],new_promB_res[i]

    new_igd_clstr = list(map(int, new_igd_clstr))

    countr = 0
    new_dist_cluster = []
    for i in range(len(new_igd_clstr)):
        if i == 0:
            countr = 1
            new_dist_cluster.append(countr)
        elif new_igd_clstr[i] - 1 != new_igd_clstr[i - 1]:
            countr = 1
            new_dist_cluster.append(countr)
        elif new_igd_clstr[i] - new_igd_clstr[i - 1] == 1:
            countr += 1
            new_dist_cluster.append(countr)
    ##    print '\n\n\nPROCESSING 4'
    ##    for i in range(20):
    ##        print new_dist_cluster[i]

    cluster_tmp = 0
    clustr_no = []
    for i in range(len(new_dist_cluster)):
        if new_dist_cluster[i] == 1:
            cluster_tmp += 1
            clustr_no.append(cluster_tmp)
        elif new_dist_cluster[i] > 1:
            cluster_tmp = cluster_tmp
            clustr_no.append(cluster_tmp)

    ##    print '\n\n\nPROCESSING 5'
    ##    for i in range(20):
    ##        print clustr_no[i]

    ##    print 'Total clusters = ', max(clustr_no), '\nLongest Cluster = :', max(new_dist_cluster),'genes long\n\n'

    if (
        len(clustr_no) == len(new_igd_clstr)
        and len(new_igd_clstr) == len(new_dist_cluster)
        and len(new_dist_cluster) == len(new_igd_scaffold_name)
        and len(new_igd_scaffold_name) == len(new_igd_geneA)
        and len(new_igd_geneA) == len(new_igd_geneB)
        and len(new_igd_geneB) == len(new_intergenic_dist)
        and len(new_intergenic_dist) == len(new_igd_stA)
        and len(new_igd_stA) == len(new_igd_stB)
        and len(new_igd_stB) == len(new_promA_res)
        and len(new_promA_res) == len(new_promB_res)
    ):
        pass
        print("Prom_IGD_Clustering___Saving...")
        zoo = zip(
            clustr_no,
            new_dist_cluster,
            new_igd_scaffold_name,
            new_igd_geneA,
            new_igd_geneB,
            new_intergenic_dist,
            new_igd_stA,
            new_igd_stB,
            new_promA_res,
            new_promB_res,
        )
        with open(output_dir + "/Operon_File.tab", "w+") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        outfile.close()
        with open(output_dir + "/IGD_PROM_CLUSTER_HDR.tab", "w+") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerow(["#", "Len", "Scf_name", "Gene1", "Gene2", "IGD", "St1", "St2", "Prom_G1", "Prom_G2"])
            writer.writerows(zoo)
        print("Done")
        outfile.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()
    print("Prom_IGD_Clustering completed")
    return


##--------------------------------------------------------------------


def Prom_clustering():
    print("Prom_clustering start")
    igd_file = output_dir + "/Gene_file_UDlen.tab"

    prom_scf_gene, prom_res = [], []
    counter = 0
    promoter_parsed_file = open(output_dir + "/Promoter_1.tab", "r")
    print("Prom_clustering___File Parsing...")
    for i in promoter_parsed_file:
        i = i.rstrip()
        counter += 1
        tmp1 = re.split(" |\):\t", i)
        prom_scf_gene.append(tmp1[0])
        prom_res.append(tmp1[-1])

    if len(prom_scf_gene) == len(prom_res):
        print("Saving...")
        zoo = zip(prom_scf_gene, prom_res)
        with open(output_dir + "/Promoter_edited.tab", "w+") as out_file:
            writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        out_file.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()

    igd_scf_nm = [x[0] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_gene = [x[1] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]
    igd_st = [x[2] for x in csv.reader(open(igd_file, "r"), delimiter="\t")]

    igd_scf_gene = []
    for i in range(len(igd_scf_nm)):
        tmp1 = igd_scf_nm[i] + "," + igd_gene[i]
        ##    tmp2 = igd_scf_nm[i]+','+igd_gene2[i]
        igd_scf_gene.append(tmp1)

    print("Prom_clustering___Promter search...")
    princ_sc_nm, res_prom, scnm_igd, gene_igd, st_igd, IGD_igd = [], [], [], [], [], []
    no_ans = "NA"
    counter = 0
    counter2 = 0
    counter3 = 0
    counter4 = 0
    counter5 = 0

    for i in igd_scf_gene:
        if i in prom_scf_gene:
            indx_igd = igd_scf_gene.index(i)
            indx_prom = prom_scf_gene.index(i)
            counter += 1
            princ_sc_nm.append(i), scnm_igd.append(igd_scf_nm[indx_igd]), gene_igd.append(
                igd_gene[indx_igd]
            ), st_igd.append(igd_st[indx_igd]), res_prom.append(prom_res[indx_prom])
        ##        print counter, i, '\t',  igd_scf_nm[indx_igd], igd_gene[indx_igd], igd_st[indx_igd], prom_res[indx_prom]
        elif i not in prom_scf_gene:
            indx_igd = igd_scf_gene.index(i)
            counter += 1
            princ_sc_nm.append(i), scnm_igd.append(igd_scf_nm[indx_igd]), gene_igd.append(
                igd_gene[indx_igd]
            ), st_igd.append(igd_st[indx_igd]), res_prom.append(no_ans)
    ##        print counter, i, '\t',  igd_scf_nm[indx_igd], igd_gene[indx_igd], igd_st[indx_igd], no_ans

    scaffold_name = scnm_igd
    geneA = gene_igd
    stA = st_igd
    prom_res = res_prom

    ##print 'scaffold_name', len(scaffold_name)
    ##print 'geneA', len(geneA)
    ##print 'stA', len(stA)
    ##print 'prom_res', len(prom_res)

    New_scaffold_name, new_geneA, new_stA, new_prom_res = [], [], [], []
    countera = 0
    counter11 = 0
    if len(scaffold_name) == len(geneA) and len(geneA) == len(stA) and len(stA) == len(prom_res):
        for i in range(len(scaffold_name)):
            if i < len(scaffold_name) - 1:

                if (
                    scaffold_name[i] == scaffold_name[i + 1]
                    and stA[i] == stA[i + 1]
                    and prom_res[i] == "Yes"
                    and prom_res[i + 1] != "Yes"
                ):
                    New_scaffold_name.append(scaffold_name[i]), new_geneA.append(geneA[i]), new_stA.append(
                        stA[i]
                    ), new_prom_res.append(prom_res[i])
                    New_scaffold_name.append(scaffold_name[i + 1]), new_geneA.append(geneA[i + 1]), new_stA.append(
                        stA[i + 1]
                    ), new_prom_res.append(prom_res[i + 1])
                ##                print scaffold_name[i], geneA[i], stA[i], prom_res[i]
                ##                print scaffold_name[i+1], geneA[i+1], stA[i], prom_res[i+1]

                elif (
                    scaffold_name[i] == scaffold_name[i + 1]
                    and stA[i] == stA[i + 1]
                    and prom_res[i] != "Yes"
                    and prom_res[i + 1] != "Yes"
                ):
                    New_scaffold_name.append(scaffold_name[i]), new_geneA.append(geneA[i]), new_stA.append(
                        stA[i]
                    ), new_prom_res.append(prom_res[i])
                    New_scaffold_name.append(scaffold_name[i + 1]), new_geneA.append(geneA[i + 1]), new_stA.append(
                        stA[i + 1]
                    ), new_prom_res.append(prom_res[i + 1])
                ##                print scaffold_name[i], geneA[i], stA[i], prom_res[i]
                ##                print scaffold_name[i+1], geneA[i+1], stA[i], prom_res[i+1]

                ##            elif scaffold_name[i] == scaffold_name[i+1] and stA[i] == stA[i+1] and prom_res[i] != 'Yes'  and prom_res[i+1] == 'Yes':
                ##                New_scaffold_name.append(scaffold_name[i]), new_geneA.append(geneA[i]), new_stA.append(stA[i]), new_prom_res.append(prom_res[i])
                ##                New_scaffold_name.append(scaffold_name[i+1]), new_geneA.append(geneA[i+1]), new_stA.append(stA[i+1]), new_prom_res.append(prom_res[i+1])
                ##                print scaffold_name[i], geneA[i], stA[i], prom_res[i]
                ##                print scaffold_name[i+1], geneA[i+1], stA[i], prom_res[i+1]
                else:
                    countera += 1
                    pass
    else:
        print("len of output is not equal")

    ##    print 'scaffold_name', len(scaffold_name)
    ##    print 'geneA', len(geneA)
    ##    print 'stA', len(stA)
    ##    print 'prom_res', len(prom_res)
    ##    print 'total entries processed ', len(prom_res)
    ##
    ##    print '\nNew_scaffold_name', len(New_scaffold_name)
    ##    print 'new_geneA', len(new_geneA)
    ##    print 'new_stA', len(new_stA)
    ##    print 'new_prom_res', len(new_prom_res)

    New_scaffold_name_uni, new_geneA_uni, new_stA_uni, new_prom_res_uni = [], [], [], []
    for i in range(len(new_geneA)):
        if i < len(new_geneA) - 1:
            if new_geneA[i] != new_geneA[i + 1]:
                new_geneA_uni.append(new_geneA[i]), New_scaffold_name_uni.append(
                    New_scaffold_name[i]
                ), new_stA_uni.append(new_stA[i]), new_prom_res_uni.append(new_prom_res[i])
            ##            print i, New_scaffold_name[i], new_geneA[i], new_stA[i], new_prom_res[i]
            else:
                pass
    ##print new_geneA_uni[:20]
    ##print '\nNew_scaffold_name_uni', len(New_scaffold_name_uni)
    ##print 'new_geneA_uni', len(new_geneA_uni)
    ##print 'new_stA_uni', len(new_stA_uni)
    ##print 'new_prom_res_uni', len(new_prom_res_uni)

    ##for i in range(len(New_scaffold_name_uni)):
    ##    print New_scaffold_name_uni[i], new_geneA_uni[i], new_stA_uni[i], new_prom_res_uni[i]

    gene_tmp = []
    for i in new_geneA_uni:
        i = i.rstrip()
        counter += 1
        tmp1 = re.split("[a-z]*_", i)
        gene_tmp.append(tmp1[2])
    ##    print gene_tmp[:20]
    gene_tmp = list(map(int, gene_tmp))

    ##print 'new_geneA_uni', len(new_geneA_uni)
    ##print 'gene_tmp', len(gene_tmp)

    cluster_len = []
    tm = 0
    for i in range(len(new_geneA_uni)):
        if (
            i == 0
            or gene_tmp[i] - 1 != gene_tmp[i - 1]
            or new_stA_uni[i] != new_stA_uni[i - 1]
            or new_prom_res_uni[i] == "Yes"
        ):
            tm = 1
            cluster_len.append(tm)
        ##        print tm, New_scaffold_name_uni[i], new_geneA_uni[i], new_stA_uni[i], new_prom_res_uni[i]
        else:
            tm += 1
            cluster_len.append(tm)
    ##        print tm, New_scaffold_name_uni[i], new_geneA_uni[i], new_stA_uni[i], new_prom_res_uni[i]

    cluster_len = list(map(int, cluster_len))
    cluster_no = []
    mp = 0
    for i in range(len(cluster_len)):
        if cluster_len[i] == 1:
            mp += 1
            cluster_no.append(mp)
        else:
            mp = mp
            cluster_no.append(mp)
    if (
        len(cluster_len) == len(cluster_no)
        and len(cluster_no) == len(New_scaffold_name_uni)
        and len(New_scaffold_name_uni) == len(new_geneA_uni)
        and len(new_geneA_uni) == len(new_stA_uni)
        and len(new_stA_uni) == len(new_prom_res_uni)
    ):
        pass
    ##        print 'print all lens equal'
    ##    for i in range(len(New_scaffold_name_uni)):
    ##        print cluster_no[i], cluster_len[i], New_scaffold_name_uni[i], new_geneA_uni[i], new_stA_uni[i], new_prom_res_uni[i]
    else:
        print("len not equal")
    ##for i in range(len(New_scaffold_name_uni)):
    ##    print cluster_no[i], cluster_len[i], New_scaffold_name_uni[i], new_geneA_uni[i], new_stA_uni[i], new_prom_res_uni[i]

    ##    print '\nNew_scaffold_name_uni', len(New_scaffold_name_uni)
    ##    print 'new_geneA_uni', len(new_geneA_uni)
    ##    print 'new_stA_uni', len(new_stA_uni)
    ##    print 'new_prom_res_uni', len(new_prom_res_uni)
    ##    print 'gene_tmp', len(gene_tmp)
    ##    print  'cluster_len', len(cluster_len)
    ##    print 'Total Entries processed: ', len(cluster_no)
    ##    print 'total clusters : ', max(cluster_no)
    ##    print 'Longest cluster : ', max(cluster_len)

    ##
    if (
        len(New_scaffold_name_uni) == len(new_geneA_uni)
        and len(new_geneA_uni) == len(new_stA_uni)
        and len(new_stA_uni) == len(new_prom_res_uni)
        and len(new_prom_res_uni) == len(cluster_len)
        and len(cluster_len) == len(cluster_no)
    ):
        pass
        print("Prom_clustering___Saving...")
        zoo = zip(cluster_no, cluster_len, New_scaffold_name_uni, new_geneA_uni, new_stA_uni, new_prom_res_uni)
        with open(output_dir + "/PROM_CLUSTER_noHDR.tab", "w+") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerows(zoo)
        outfile.close()
    else:
        print("Number of entries for each gene are not equal, please correct the data and run again")
        exit()
    print("Prom_clustering completed")

    return


def config_file_check():
    global NNPP2_path
    print("config_file_check start")
    My_prog_path = os.path.dirname(os.path.abspath(__file__))
    print(My_prog_path)
    if os.path.exists(My_prog_path + "/NNPP2.2/bin/fa2TDNNpred-PRO-static.linux") == True:
        NNPP2_path = My_prog_path + "/NNPP2.2"
    else:
        print("please check the path for NNPP2.2")
        sys.exit(1)

    return


##def Processing_details(file_name, start_time, total_loops):
def Processing_details(file_name, start_time, counter):
    ##    print "Total entries processed", len(file_name)
    ##    print 'Total Promoter prediction attempted', counter
    end_time = time.time()
    elapsed = end_time - start_time
    print("Total time elapsed", elapsed)
    print("Please check ", output_dir)

    try:
        os.remove(output_dir + "/Gene_coord.tab")
        os.remove(output_dir + "/DIR_Proximons_hdr.tab")
        os.remove(output_dir + "/DIR_Proximons_temp.tab")
        os.remove(output_dir + "/Downstream_seq.fasta")
        os.remove(output_dir + "/Upstream_seq.fasta")
        os.remove(output_dir + "/Gene_ups_Coordinates.tab")
        os.remove(output_dir + "/Gene_dss_Coordinates.tab")
        os.remove(output_dir + "/Promoter_1.tab")
        os.remove(output_dir + "/IGD_PROM_RES.tab")
        os.remove(output_dir + "/PROM_CLUSTER_noHDR.tab")
        os.remove(output_dir + "/Promoter_1.tab")
        os.remove(output_dir + "/promoter_prediction.txt")
        os.remove(output_dir + "/IGD_PROM_CLUSTER_HDR.tab")
        os.remove(output_dir + "/promoter_prediction.txt")
        os.remove(output_dir + "/Gene_file_detailed.tab")
        os.remove(output_dir + "/Gene_file_edited.tab")
    except OSError:
        pass


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("I got interrupted. :-( Bye!\n")
        sys.exit(0)
