import os
import re
import logging
import shutil
from datetime import datetime
import multiprocessing
from functools import partial

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import pan_genome.utils as utils


def create_pro_file_for_each_cluster(samples, gene_to_cluster, collection_dir):
    """
    Create fasta files, containing protein sequences for each gene cluster.

    Parameters
    ----------
    samples : list of dict
        list of samples information {id: , gff_file: , assembly: }
    gene_to_cluster : dict
        find the cluster of a specific gene
        {gene_id : cluster_id}
    collection_dir
        the directory of collection
    """
    starttime = datetime.now()
    clusters_dir = os.path.join(collection_dir, "clusters")
    # remove old sequence file
    os.system(f'find {clusters_dir} -name "*.faa" -delete')

    for sample in samples:
        sample_id = sample["id"]
        faa_file = os.path.join(collection_dir, "samples", sample_id, sample_id + ".faa")
        with open(faa_file, "r") as in_fh:
            last_seq_id = None
            seq_lines = []
            for line in in_fh:
                line = line.rstrip()
                result = re.match(r"^>(\S+)", line)
                if result != None:
                    seq_id = result.group(1)
                    if last_seq_id != None:
                        cluster_id = str(gene_to_cluster[last_seq_id])
                        cluster_dir = os.path.join(collection_dir, "clusters", cluster_id)
                        if not os.path.exists(cluster_dir):
                            os.mkdir(cluster_dir)
                        cluster_file = os.path.join(cluster_dir, cluster_id + ".faa")
                        with open(cluster_file, "a") as out_fh:
                            out_fh.write(">" + last_seq_id + "\n")
                            out_fh.write("\n".join(seq_lines) + "\n")
                    last_seq_id = seq_id
                    seq_lines = []
                else:
                    seq_lines.append(line)

            cluster_id = str(gene_to_cluster[last_seq_id])
            cluster_dir = os.path.join(collection_dir, "clusters", cluster_id)
            if not os.path.exists(cluster_dir):
                os.mkdir(cluster_dir)
            cluster_file = os.path.join(cluster_dir, cluster_id + ".faa")
            with open(cluster_file, "a") as out_fh:
                out_fh.write(">" + last_seq_id + "\n")
                out_fh.write("\n".join(seq_lines) + "\n")

    elapsed = datetime.now() - starttime
    logging.info("Create protein sequence file for each gene cluster " f"-- time taken {str(elapsed)}")


def create_poa(cluster_id, collection_dir, baseDir):
    """
    Run abPOA to create mutiple sequence alignment of a cluster.
    Call the consensus sequence (the representative sequence).

    Parameters
    ----------
    cluster_id : int
        the id of cluster
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    """
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, "clusters", cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + ".faa")

    matrix_file = os.path.join(baseDir, "BLOSUM62.mtx")
    result_file = os.path.join(cluster_dir, cluster_id + ".result")
    cmd = f"abpoa {seq_file} -o {result_file} -r2 -t {matrix_file}" " -O 11,0 -E 1,0 -p -c 2> /dev/null"
    if not os.path.exists(result_file):
        utils.run_command(cmd)


def create_poa_in_parallel(clusters_id_list, collection_dir, baseDir, threads):
    """
    Run abPOA in parallel.

    Parameters
    ----------
    clusters_id_list : list
        list of cluster id
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    out_dir : path
        output directory
    threads : int
        number of threads
    """
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(partial(create_poa, collection_dir=collection_dir, baseDir=baseDir), clusters_id_list)

    elapsed = datetime.now() - starttime
    logging.info("Create multiple sequence alignment by abPOA " f"-- time taken {str(elapsed)}")


def add_poa(cluster_id, collection_dir, baseDir):
    """
    Use abPOA to add new sequences in to an existing MSA of the previous
    cluster. Then, get the consensus sequence (the representative sequence).

    Parameters
    ----------
    cluster_id : int
        the id of cluster
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    """
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, "clusters", cluster_id)
    seq_file = os.path.join(cluster_dir, cluster_id + ".faa")
    msa_file = os.path.join(cluster_dir, cluster_id + ".msa")
    matrix_file = os.path.join(baseDir, "BLOSUM62.mtx")
    result_file = os.path.join(cluster_dir, cluster_id + ".result")
    # if there is no new sequences, not run abPOA again. Parse the
    # previous result file to get the consensus sequence.
    cmd = (
        f"abpoa {seq_file} -i {msa_file} -o {result_file} -r2 "
        f"-t {matrix_file} -O 11,0 -E 1,0 -p -c -m 1 2> /dev/null"
    )
    if os.path.isfile(seq_file) and not os.path.exists(result_file):
        utils.run_command(cmd)


def add_poa_in_parallel(clusters_id_list, collection_dir, baseDir, threads):
    """
    Run abPOA in parallel to add new sequences into previous clusters.

    Parameters
    ----------
    clusters_id_list : list
        list of cluster id
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    threads : int
        number of threads
    """
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        pool.map(partial(add_poa, collection_dir=collection_dir, baseDir=baseDir), clusters_id_list)

    elapsed = datetime.now() - starttime
    logging.info("Add new sequences to previous multiple sequence alignment " f"by abPOA -- time taken {str(elapsed)}")


def extract_consensus_seq(cluster_id, collection_dir):
    """
    Extract the consensus sequence (the representative sequence).

    Parameters
    ----------
    cluster_id : int
        the id of cluster
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix

    Returns
    -------
    cluster_id : int
        the id of this cluster
    cons_seq : str
        consensus sequence of this cluster
    """
    cluster_id = str(cluster_id)
    cluster_dir = os.path.join(collection_dir, "clusters", cluster_id)
    result_file = os.path.join(cluster_dir, cluster_id + ".result")
    seq_file = os.path.join(cluster_dir, cluster_id + ".faa")
    cons_seq = ""
    msa_file = os.path.join(cluster_dir, cluster_id + ".msa")
    if not os.path.exists(result_file):
        # if there is no new sequences, the is no result file
        return cluster_id, None

    with open(result_file, "r") as in_fh, open(msa_file, "w") as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        for seq_record in SeqIO.parse(in_fh, "fasta"):
            if seq_record.id == "Consensus_sequence":
                cons_seq = str(seq_record.seq)
            else:
                fasta_out.write_record(seq_record)
    os.remove(result_file)
    os.remove(seq_file)

    return cluster_id, cons_seq


def combine_representative_seqs(clusters_id_list, collection_dir, threads):
    """
    Combine all consensus sequences of gene clusters into a fasta file, named
    reference_pangenome.fasta.

    Parameters
    ----------
    clusters_id_list : list
        list of cluster id
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    threads : int
        number of threads

    Returns
    -------
    representative_fasta : path
        a fasta file contain all representative sequences of clusters.
    """
    starttime = datetime.now()

    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(partial(extract_consensus_seq, collection_dir=collection_dir), clusters_id_list)

    representative_fasta = os.path.join(collection_dir, "reference_pangenome.fasta")
    if os.path.exists(representative_fasta):
        previous_representative = {}
        with open(representative_fasta, "r") as in_fh:
            for seq_record in SeqIO.parse(in_fh, "fasta"):
                cluster_id = seq_record.id
                seq = str(seq_record.seq)
                previous_representative[cluster_id] = seq

    temp_representative_fasta = os.path.join(collection_dir, "reference_pangenome.fasta.temp")
    with open(temp_representative_fasta, "w") as out_fh:
        for cluster_id, cons_seq in results:
            # if there is no new sequences, the consensus will not change
            if cons_seq == None:
                cons_seq = previous_representative[cluster_id]
            # remove dash character in consensus sequences
            cons_seq = re.sub("-", "", cons_seq)
            new_record = SeqRecord(Seq(cons_seq), id=cluster_id, description="")
            SeqIO.write(new_record, out_fh, "fasta")

    shutil.move(temp_representative_fasta, representative_fasta)

    elapsed = datetime.now() - starttime
    logging.info("Create representative fasta" f"-- time taken {str(elapsed)}")

    return representative_fasta


def create_msa_init_pipeline(clusters, samples, collection_dir, baseDir, threads):
    """
    Create multiple sequence alignment for Init pipeline.

    Parameters
    ----------
    clusters : list of list
        list of sequence IDs of each cluster
    samples : list of dict
        list of samples information {id: , gff_file: , assembly: }
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    threads : int
        number of threads

    Returns
    -------
    representative_fasta : path
        a fasta file contain all representative sequences of clusters.
    """
    # skip this step if it has been done.
    finished = os.path.join(collection_dir, "temp", "alignment.done")
    if os.path.exists(finished):
        representative_fasta = os.path.join(collection_dir, "reference_pangenome.fasta")
        logging.info("Skip alignment")
        return representative_fasta

    clusters_dir = os.path.join(collection_dir, "clusters")
    if not os.path.exists(clusters_dir):
        os.mkdir(clusters_dir)
    else:
        shutil.rmtree(clusters_dir)
        os.mkdir(clusters_dir)
    # create dictionary to find the cluster of a specific gene
    gene_to_cluster_id = {gene: i for i, cluster in enumerate(clusters) for gene in cluster}
    clusters_id_list = range(0, len(clusters))

    # create protein sequences for each clusters
    create_pro_file_for_each_cluster(samples, gene_to_cluster_id, collection_dir)

    # create msa
    create_poa_in_parallel(clusters_id_list, collection_dir, baseDir, threads)

    representative_fasta = combine_representative_seqs(clusters_id_list, collection_dir, threads)

    # create a file to indicate that this step is done.
    os.system(f"touch {finished}")

    return representative_fasta


def create_msa_add_pipeline(previous_clusters, new_clusters, new_samples, collection_dir, baseDir, threads):
    """
    Create multiple sequence alignment for Add pipeline.
    + add new sequences to previous MSA
    + create MSA for new clusters

    Parameters
    ----------
    previous_clusters : list of list
        list of sequence IDs of previous clusters
    new_clusters : list of list
        list of sequence IDs of new clusters
    new_samples : list of dict
        list of new samples information {id: , gff_file: , assembly: }
    collection_dir : path
        collection directory
    baseDir : path
        panta directory, contain BLOSUM62 matrix
    threads : int
        number of threads

    Returns
    -------
    representative_fasta : path
        a fasta file contain all representative sequences of clusters.
    """
    # skip this step if it has been done.
    finished = os.path.join(collection_dir, "temp", "alignment.done")
    if os.path.exists(finished):
        representative_fasta = os.path.join(collection_dir, "reference_pangenome.fasta")
        logging.info("Skip alignment")
        return representative_fasta
    # find the cluster id for a specific gene
    gene_to_cluster_id = {}
    for i, cluster in enumerate(previous_clusters):
        for gene in cluster:
            gene_to_cluster_id[gene] = i
    # continue the id of old clusters
    for i, cluster in enumerate(new_clusters, len(previous_clusters)):
        for gene in cluster:
            gene_to_cluster_id[gene] = i

    # get new gene sequences for each clusters
    create_pro_file_for_each_cluster(new_samples, gene_to_cluster_id, collection_dir)

    # add new gene to existing clusters
    cluster_id_list = range(0, len(previous_clusters))
    add_poa_in_parallel(cluster_id_list, collection_dir, baseDir, threads)

    # create msa for the new clusters
    cluster_id_list = range(len(previous_clusters), len(previous_clusters) + len(new_clusters))
    create_poa_in_parallel(cluster_id_list, collection_dir, baseDir, threads)
    cluster_id_list = range(0, len(previous_clusters) + len(new_clusters))

    representative_fasta = combine_representative_seqs(cluster_id_list, collection_dir, threads)

    # create a file to indicate that this step is done.
    os.system(f"touch {finished}")

    return representative_fasta
