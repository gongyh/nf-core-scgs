import os
import re
import logging
import multiprocessing
from functools import partial
from datetime import datetime
import gzip
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pan_genome.utils as utils
import pan_genome.output as output

logger = logging.getLogger(__name__)


def parse_gff_file(ggf_file, sample_dir, sample_id, add_strand_2gene):
    """
    Parse gff file.
    - Filter out too short gene (less than 120 nucleotides).
    - Extract gene information, put into gene dictionary.
    - Write a bed file, containing location of each gene.
    - Write a fasta file, which is genome assembly.

    Parameters
    ----------
    gff_file : path
        path to gff file
    sample_dir : path
        directory of output
    sample_id : str
        sample ID

    Returns
    -------
    bed_file : path
        location of genes in BED format
    assembly_file : path
        genome assembly in FASTA format
    gene dictionary : dict
        contain information of each gene of one sample
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    """
    bed_file = os.path.join(sample_dir, sample_id + ".bed")
    assembly_file = os.path.join(sample_dir, sample_id + ".fasta")
    gene_dictionary = {}
    gene_position = {}
    found_fasta = False
    suffix = 1

    if ggf_file.endswith(".gz"):
        in_fh = gzip.open(ggf_file, "rt")
    else:
        in_fh = open(ggf_file, "r")

    with open(bed_file, "w") as bed_fh, open(assembly_file, "w") as fna_fh:
        for line in in_fh:
            if found_fasta == True:
                # write genome assembly
                fna_fh.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                continue
            if re.match(r"^#", line) != None:
                continue
            line = line.rstrip("\n")
            cells = line.split("\t")
            if cells[2] != "CDS":
                continue

            start = int(cells[3])
            end = int(cells[4])
            length = end - start + 1
            # if length < 120: # filter out gene less 120 nu
            #     print('short')
            #     continue
            if length % 3 != 0:
                continue
            contig = cells[0]
            trand = cells[6]
            tags = cells[8].split(";")
            gene_id = None
            gene_name = ""
            gene_product = ""
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    gene_id = ID.group(1)
                    # gene_id = re.sub(r'\W', '_', gene_id)
                    continue

                gene = re.match(r"^gene=(.+)", tag)
                if gene != None:
                    gene_name = gene.group(1)
                    gene_name = re.sub(r"\W", "_", gene_name)
                    gene_name = re.sub(r"_\d+$", "", gene_name)
                    continue

                product = re.match(r"^product=(.+)", tag)
                if product != None:
                    gene_product = product.group(1)
            if gene_id == None:
                continue

            gene_id += f"_{sample_id}_{str(suffix)}"
            suffix += 1

            # add strand information to gene name (default = False)
            if add_strand_2gene:
                gene_id += "@" + str(length) + "@"
                gene_id += trand

            # create bed file
            row = [contig, str(start - 1), str(end), gene_id, "1", trand]
            bed_fh.write("\t".join(row) + "\n")
            # add to gene_dictionary
            gene_dictionary[gene_id] = (sample_id, contig, length, gene_name, gene_product)
            gene_position.setdefault(contig, []).append(gene_id)

    in_fh.close()

    return bed_file, assembly_file, gene_dictionary, gene_position


def process_single_sample_gff(sample, out_dir, table, add_strand_2gene):
    """
    Process sample if input data is genome annotation.

    Parse GFF file, extract gene sequences by BEDTools, translate
    nucleotide to protein.

    Parameters
    ----------
    sample : dict
        sample information {id: , gff_file: , assembly: }
    out_dir : path
        output directory
    table : int
        codon table

    Returns
    -------
    gene dictionary : dict
        contain information of each gene of one sample
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    """
    sample_id = sample["id"]
    sample_dir = os.path.join(out_dir, "samples", sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

    # parse gff file
    bed_file, assembly_file, gene_dictionary, gene_position = parse_gff_file(
        ggf_file=sample["gff_file"], sample_dir=sample_dir, sample_id=sample_id, add_strand_2gene=add_strand_2gene
    )
    # if input has a separated genome assembly, use it.
    if sample["assembly"] != None:
        os.remove(assembly_file)  # remove an empty file
        assembly_file = sample["assembly"]

    # if result file exists, skip
    faa_file = os.path.join(sample_dir, sample_id + ".faa")
    if os.path.isfile(faa_file):
        return gene_dictionary

    # extract nucleotide region
    fna_file = os.path.join(sample_dir, sample_id + ".fna")
    cmd = f"bedtools getfasta -s -fi {assembly_file} -bed {bed_file} " f"-fo {fna_file} -name > /dev/null 2>&1"
    utils.run_command(cmd)

    # translate nucleotide to protein
    utils.translate_protein(nu_fasta=fna_file, pro_fasta=faa_file, table=table)

    if sample["assembly"] == None:
        os.remove(assembly_file)
    os.remove(bed_file)
    os.remove(assembly_file + ".fai")
    os.remove(fna_file)

    return gene_dictionary, gene_position


def process_single_sample_fasta(sample, out_dir, table, add_strand_2gene):
    """
    Process sample if input is genome assembly.

    + Predict gene by Prodigal.
    + Filter out gene sequences.
    + Rename gene id and write a new protein sequences fasta file.
    + Extract gene information, and put into gene dictionary.

    Parameters
    ----------
    sample : dict
        sample information {id: , gff_file: , assembly: }
    out_dir : path
        output directory
    table : int
        codon table

    Returns
    -------
    gene dictionary : dict
        contain information of each gene of one sample
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    """

    sample_id = sample["id"]
    sample_dir = os.path.join(out_dir, "samples", sample_id)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

    # gene prediction
    assembly_file = sample["assembly"]
    faa_file = os.path.join(sample_dir, sample_id + ".original.faa")
    if assembly_file.endswith(".gz"):
        cmd = f"zcat {assembly_file} | " f"prodigal -a {faa_file} -g {table} -c -m -q -o /dev/null"
    else:
        cmd = f"prodigal -i {assembly_file} -a {faa_file} -g {table} " f"-c -m -q -o /dev/null"
    # if result files exists, skip Prodigal
    if not os.path.isfile(faa_file):
        utils.run_command(cmd)

    # rename gene id and extract coordinates
    gene_dictionary = {}
    gene_position = {}
    rewrite_faa_file = os.path.join(sample_dir, sample_id + ".faa")
    count = 1
    with open(faa_file, "r") as in_fh, open(rewrite_faa_file, "w") as out_fh:
        for record in SeqIO.parse(in_fh, "fasta"):
            contig = record.id.rsplit("_", 1)[0]
            cells = record.description.split(" # ")
            pro = str(record.seq)

            # filter short genes
            start = int(cells[1])
            end = int(cells[2])
            trand = cells[3]
            length = end - start + 1
            # if length < 120:
            #     # logger.info('Short gene')
            #     continue
            if length % 3 != 0:
                continue
            # filter seq with premature codon
            results = re.findall(r"\*", pro)
            if len(results) > 1:
                # logger.info('Have premature codon')
                continue

            # filter seq lacking start and stop codon
            if pro[0] != "M" and pro[-1] != "*":
                # logger.info('Lack both start and stop codon')
                continue

            # filter seq which has more than 5% of unknown
            results = re.findall(r"X", pro)
            if len(results) / len(pro) > 0.05:
                # logger.info('Too many unknowns')
                continue

            gene_id = sample_id + f"_{str(count)}"
            # add strand information to gene name (default = False)
            if add_strand_2gene:
                gene_id += "@" + str(length) + "@"
                gene_id += trand

            count += 1
            desc = "{}~~~{}~~~{}~~~{}".format(contig, start, end, trand)
            new_record = SeqRecord(record.seq, id=gene_id, description=desc)
            SeqIO.write(new_record, out_fh, "fasta")

            # add to gene_dictionary
            gene_dictionary[gene_id] = (sample_id, contig, length)
            gene_position.setdefault(contig, []).append(gene_id)

    return gene_dictionary, gene_position


def extract_proteins(samples, out_dir, args):
    """
    Process samples in parallel. Combine gene_dictionary of all samples.

    Parameters
    ----------
    samples : list of dict
        list of samples
    out_dir : path
        output directory
    args : object
        Command-line input arguments

    Returns
    -------
    gene dictionary : dict
        contain information of each gene of all samples
        {gene_id: (sample_id, contig, length, gene_name, gene_product)}
    """
    starttime = datetime.now()

    with multiprocessing.Pool(processes=args.threads) as pool:
        if args.fasta == None:
            results = pool.map(
                partial(
                    process_single_sample_gff, out_dir=out_dir, table=args.table, add_strand_2gene=args.addstrand2gene
                ),
                samples,
            )
        else:
            results = pool.map(
                partial(
                    process_single_sample_fasta, out_dir=out_dir, table=args.table, add_strand_2gene=args.addstrand2gene
                ),
                samples,
            )

    gene_dictionary = {}
    gene_position = {}
    for sample, result in zip(samples, results):
        gene_dictionary.update(result[0])
        gene_position[sample["id"]] = result[1]

    output.write_gene_position(gene_position, out_dir)

    elapsed = datetime.now() - starttime
    logging.info(f"Extract protein -- time taken {str(elapsed)}")

    return gene_dictionary, gene_position


def combine_proteins(collection_dir, out_dir, samples, timing_log, resume):
    """
    Combine protein sequences of all samples into one file.

    Parameters
    ----------
    collection_dir : path
        collection directory
    out_dir : path
        directory of output file
    samples : list of dict
        information of each sample
    timing_log : path
        path to time.log
    resume : list
        A boolean inside a list
        If True, resume previous analysis

    Returns
    -------
    path
        path of output file
    """
    # starttime = datetime.now()
    combined_faa_file = os.path.join(out_dir, "combined.faa")
    if os.path.isfile(combined_faa_file) and resume[0] == True:
        logging.info(f"Resume - Combine protein")
        return combined_faa_file
    else:
        resume[0] = False

    protein_files = os.path.join(out_dir, "protein.txt")
    with open(protein_files, "w") as fh:
        for sample in samples:
            sample_id = sample["id"]
            faa_file = os.path.join(collection_dir, "samples", sample_id, sample_id + ".faa")
            if os.path.isfile(faa_file):
                fh.write(faa_file + "\n")
            else:
                raise Exception(f"{faa_file} does not exist")

    cmd = f"cat {protein_files} | xargs cat > {combined_faa_file}"
    utils.run_command(cmd, timing_log)

    # elapsed = datetime.now() - starttime
    # logging.info(f'Combine protein -- time taken {str(elapsed)}')
    return combined_faa_file
