import os
import re
import logging
import gzip
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import pan_genome.utils as utils


logger = logging.getLogger(__name__)


def find_paralogs(cluster, gene_dictionary):
    samples = {}
    for gene_id in cluster:
        sample_id = gene_dictionary[gene_id][0]
        samples.setdefault(sample_id, []).append(gene_id)

    # pick paralogs with the smallest number of genes
    smallest_number = 1000000
    paralog_genes = None
    for sample_id in samples:
        genes = samples[sample_id]
        count = len(genes)
        if count > 1 and count < smallest_number:
            paralog_genes = genes
            smallest_number = count

    return paralog_genes


def get_neighbour_genes(gene_dictionary, gene_position):
    gene_neighbour_dict = {}
    for gene_id in gene_dictionary:
        contig = gene_dictionary[gene_id][1]
        sample_id = gene_dictionary[gene_id][0]
        genes_of_contig = gene_position[sample_id][contig]
        index = genes_of_contig.index(gene_id)
        pre_index = index - 5
        post_index = index + 6
        if pre_index < 0:
            pre_index = 0
        length_of_contig = len(genes_of_contig)
        if post_index >= length_of_contig:
            post_index = length_of_contig
        neighbour_genes = genes_of_contig[pre_index:index] + genes_of_contig[index + 1 : post_index]
        gene_neighbour_dict[gene_id] = neighbour_genes

    return gene_neighbour_dict


def create_orthologs(cluster, paralog_genes, gene_neighbour_dict, gene_to_cluster_index):
    # find cluster indices of all the neighbour genes
    # of each paralog gene
    cluster_indices_around_paralogs = []
    for p in paralog_genes:
        neighbours_of_p = gene_neighbour_dict[p]
        cluster_indices_around_p = set()
        for neighbour_gene in neighbours_of_p:
            try:
                cluster_index = gene_to_cluster_index[neighbour_gene]
                cluster_indices_around_p.add(cluster_index)
            except:
                continue
        cluster_indices_around_paralogs.append(cluster_indices_around_p)

    # create data structure to hold new clusters
    new_clusters = [[p] for p in paralog_genes]
    # extra "leftovers" list to gather genes
    # that don't share CGN with any paralog gene
    new_clusters.append([])

    # add other members of the cluster to their closest match
    for g in cluster:
        if g in paralog_genes:
            continue

        neighbour_genes_of_g = gene_neighbour_dict[g]
        if len(neighbour_genes_of_g) == 0:
            new_clusters[-1].append(g)
            continue

        # find paralog gene which is closest match with g
        best_score = 0
        best_score_index = -1  # -1 is the index of "leftovers" list
        for p, _ in enumerate(paralog_genes):
            cluster_indices_around_p = cluster_indices_around_paralogs[p]
            score_of_p = 0
            for neighbour_gene in neighbour_genes_of_g:
                try:
                    cluster_index = gene_to_cluster_index[neighbour_gene]
                except:
                    continue
                if cluster_index in cluster_indices_around_p:
                    score_of_p += 1
            score_of_p = score_of_p / len(neighbour_genes_of_g)
            if score_of_p > best_score:
                best_score = score_of_p
                best_score_index = p

        new_clusters[best_score_index].append(g)

    # check for "leftovers", remove if absent
    if len(new_clusters[-1]) == 0:
        del new_clusters[-1]

    return new_clusters


def split_paralogs(gene_dictionary, gene_position, unsplit_clusters, split):
    if split == False:
        return unsplit_clusters

    starttime = datetime.now()

    gene_neighbour_dict = get_neighbour_genes(gene_dictionary, gene_position)

    clusters_not_paralogs = set()
    # run iteratively
    out_clusters = unsplit_clusters
    for i in range(100000):
        stime = datetime.now()
        in_clusters = out_clusters
        out_clusters = []
        any_paralogs = 0
        # convert in_clusters so we can find the cluster index of gene
        gene_to_cluster_index = {gene: index for index, genes in enumerate(in_clusters) for gene in genes}

        for cluster in in_clusters:
            if len(cluster) == 1:
                out_clusters.append(cluster)
                continue
            first_gene = cluster[0]
            if first_gene in clusters_not_paralogs:
                out_clusters.append(cluster)
                continue

            # check paralogs
            paralog_genes = find_paralogs(cluster, gene_dictionary)

            if paralog_genes == None:
                clusters_not_paralogs.add(first_gene)
                out_clusters.append(cluster)
                continue

            # split paralogs
            orthologs_clusters = create_orthologs(cluster, paralog_genes, gene_neighbour_dict, gene_to_cluster_index)
            out_clusters.extend(orthologs_clusters)
            any_paralogs = 1

        # check if next iteration is required
        if any_paralogs == 0:
            break
        elapsed = datetime.now() - stime
        logging.info(f"Split paralogs iterate {i}-- time taken {str(elapsed)}")
    split_clusters = out_clusters

    elapsed = datetime.now() - starttime
    logging.info(f"Split paralogs -- time taken {str(elapsed)}")
    return split_clusters


def create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, out_dir):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, "clusters")
    pan_ref_file = os.path.join(out_dir, "pan_genome_reference.fna")
    with open(pan_ref_file, "w") as ref_fh:
        for sample in samples:
            sample_id = sample["id"]
            fna_file = os.path.join(out_dir, "samples", sample_id, sample_id + ".fna")
            with open(fna_file, "r") as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if re.match(r"^>", line) != None:
                        line = re.sub(r"\([-+]\)", "", line)
                        result = re.match(r"^>([^:]+)", line)
                        seq_id = result.group(1)
                    else:
                        ls = [line[i : i + 60] for i in range(0, len(line), 60)]
                        if seq_id in gene_to_cluster_name:
                            cluster_name = gene_to_cluster_name[seq_id]
                            cluster_file = os.path.join(clusters_dir, cluster_name, cluster_name + ".fna")
                            with open(cluster_file, "a") as out_fh:
                                out_fh.write(">" + seq_id + "\n")
                                out_fh.write("\n".join(ls) + "\n")
                        if seq_id in pan_ref_list:
                            ref_fh.write(">" + seq_id + "\n")
                            ref_fh.write("\n".join(ls) + "\n")

    elapsed = datetime.now() - starttime
    logging.info("Create nucleotide sequence file for each gene cluster " f"-- time taken {str(elapsed)}")


def create_pro_file_for_each_cluster(samples, gene_to_cluster_name, out_dir):
    starttime = datetime.now()

    for sample in samples:
        sample_id = sample["id"]
        faa_file = os.path.join(out_dir, "samples", sample_id, sample_id + ".faa")
        with open(faa_file, "r") as in_fh:
            last_seq_id = None
            line_list = []
            for line in in_fh:
                line = line.rstrip()
                result = re.match(r"^>(.+)", line)
                if result != None:
                    seq_id = result.group(1)
                    if last_seq_id != None:
                        cluster_name = gene_to_cluster_name[last_seq_id]
                        cluster_file = os.path.join(out_dir, "clusters", cluster_name, cluster_name + ".faa")
                        with open(cluster_file, "a") as out_fh:
                            out_fh.write(">" + last_seq_id + "\n")
                            out_fh.write("\n".join(line_list) + "\n")
                    last_seq_id = seq_id
                    line_list = []
                else:
                    line_list.append(line)

            cluster_name = gene_to_cluster_name[last_seq_id]
            cluster_file = os.path.join(out_dir, "clusters", cluster_name, cluster_name + ".faa")
            with open(cluster_file, "a") as out_fh:
                out_fh.write(">" + last_seq_id + "\n")
                out_fh.write("\n".join(line_list) + "\n")

    elapsed = datetime.now() - starttime
    logging.info(f"Create protein sequence file for each gene cluster " f"-- time taken {str(elapsed)}")


def run_mafft_protein_alignment(annotated_clusters, out_dir, threads, timing_log):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, "clusters")

    cmds_file = os.path.join(clusters_dir, "pro_align_cmds")
    with open(cmds_file, "w") as cmds:
        for cluster_name in annotated_clusters:
            cluster_dir = os.path.join(clusters_dir, cluster_name)
            gene_aln_file = os.path.join(cluster_dir, cluster_name + ".faa.aln.gz")
            gene_seq_file = os.path.join(cluster_dir, cluster_name + ".faa")
            if not os.path.isfile(gene_seq_file):
                logger.info("{} does not exist".format(gene_seq_file))
                continue
            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} " f"| gzip > {gene_aln_file}"
            cmds.write(cmd + "\n")

    cmd = f"parallel --progress -j {threads} -a {cmds_file}"
    utils.run_command(cmd, timing_log)

    elapsed = datetime.now() - starttime
    logging.info(f"Run protein alignment -- time taken {str(elapsed)}")


def run_mafft_nucleotide_alignment(annotated_clusters, out_dir, threads, timing_log):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, "clusters")

    cmds_file = os.path.join(clusters_dir, "nu_align_cmds")
    with open(cmds_file, "w") as cmds:
        for cluster_name in annotated_clusters:
            cluster_dir = os.path.join(clusters_dir, cluster_name)
            gene_aln_file = os.path.join(cluster_dir, cluster_name + ".fna.aln.gz")
            gene_seq_file = os.path.join(cluster_dir, cluster_name + ".fna")
            if not os.path.isfile(gene_seq_file):
                logger.info("{} does not exist".format(gene_aln_file))
                continue
            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} " f"| gzip > {gene_aln_file} && rm {gene_seq_file}"
            cmds.write(cmd + "\n")

    cmd = f"parallel --progress -j {threads} -a {cmds_file}"
    utils.run_command(cmd, timing_log)

    elapsed = datetime.now() - starttime
    logging.info(f"Run nucleotide alignment -- time taken {str(elapsed)}")


def create_nucleotide_alignment(annotated_clusters, out_dir):
    starttime = datetime.now()

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(out_dir, "clusters", cluster_name)

        protein_aln_file = os.path.join(cluster_dir, cluster_name + ".faa.aln.gz")
        if not os.path.isfile(protein_aln_file):
            logger.info("{} does not exist".format(protein_aln_file))
            continue
        protein_dict = {}
        with gzip.open(protein_aln_file, "rt") as fh:
            for seq_record in SeqIO.parse(fh, "fasta"):
                protein_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_seq_file = os.path.join(cluster_dir, cluster_name + ".fna")
        if not os.path.isfile(nucleotide_seq_file):
            logger.info("{} does not exist".format(nucleotide_seq_file))
            continue
        nucleotide_dict = {}
        for seq_record in SeqIO.parse(nucleotide_seq_file, "fasta"):
            nucleotide_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_aln_file = os.path.join(cluster_dir, cluster_name + ".fna.aln.gz")
        with gzip.open(nucleotide_aln_file, "wt") as fh:
            for seq_id in protein_dict.keys():
                protein = protein_dict[seq_id]
                nucleotide = nucleotide_dict[seq_id]
                result = ""
                codon_pos = 0
                for c in protein:
                    if c == "-":
                        result += "---"
                    else:
                        result += nucleotide[codon_pos * 3 : codon_pos * 3 + 3]
                        codon_pos += 1
                new_record = SeqRecord(Seq(result), id=seq_id, description="")
                SeqIO.write(new_record, fh, "fasta")

        os.remove(nucleotide_seq_file)
        os.remove(os.path.join(cluster_dir, cluster_name + ".faa"))

    elapsed = datetime.now() - starttime
    logging.info(f"Create  nucleotide alignment -- time taken {str(elapsed)}")


def create_core_gene_alignment(annotated_clusters, gene_dictionary, samples, out_dir):
    starttime = datetime.now()

    clusters_dir = os.path.join(out_dir, "clusters")
    seq_dict = {}
    for sample in samples:
        sample_id = sample["id"]
        seq_dict[sample_id] = ""

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(clusters_dir, cluster_name)

        sample_list = set()
        skip = False
        for gene in annotated_clusters[cluster_name]["gene_id"]:
            sample_id = gene_dictionary[gene][0]
            if sample_id not in sample_list:
                sample_list.add(sample_id)
            else:
                skip = True
        if len(sample_list) < len(samples):
            skip = True

        if skip == True:
            continue

        nucleotide_aln_file = os.path.join(cluster_dir, cluster_name + ".fna.aln.gz")
        if not os.path.isfile(nucleotide_aln_file):
            logger.info("{} does not exist".format(nucleotide_aln_file))
            continue

        cluster_dict = {}
        with gzip.open(nucleotide_aln_file, "rt") as fh:
            for seq_record in SeqIO.parse(fh, "fasta"):
                sample_name = gene_dictionary[seq_record.id][0]
                cluster_dict[sample_name] = str(seq_record.seq)

        for sample_name in cluster_dict:
            seq_dict[sample_name] += cluster_dict[sample_name]

    core_gene_aln_file = os.path.join(out_dir, "core_gene_alignment.aln.gz")
    with gzip.open(core_gene_aln_file, "wt") as fh:
        for sample in seq_dict:
            new_record = SeqRecord(Seq(seq_dict[sample]), id=sample, description="")
            SeqIO.write(new_record, fh, "fasta")

    elapsed = datetime.now() - starttime
    logging.info(f"Create core gene alignment -- time taken {str(elapsed)}")


def run_gene_alignment(annotated_clusters, gene_dictionary, samples, collection_dir, alignment, threads, timing_log):
    if alignment == None:
        return

    gene_to_cluster_name = {}
    pan_ref_list = set()

    clusters_dir = os.path.join(collection_dir, "clusters")
    if os.path.exists(clusters_dir):
        shutil.rmtree(clusters_dir)
        os.mkdir(clusters_dir)
    else:
        os.mkdir(clusters_dir)

    for cluster_name in annotated_clusters:
        cluster_dir = os.path.join(collection_dir, "clusters", cluster_name)
        if not os.path.exists(cluster_dir):
            os.mkdir(cluster_dir)
        length_max = 0
        representative = None
        for gene in annotated_clusters[cluster_name]["gene_id"]:
            gene_to_cluster_name[gene] = cluster_name
            length = gene_dictionary[gene][2]
            if length > length_max:
                representative = gene
                length_max = length
        pan_ref_list.add(representative)

    if alignment == "protein":
        pa.create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir)
        pa.create_pro_file_for_each_cluster(samples, gene_to_cluster_name, collection_dir)
        pa.run_mafft_protein_alignment(annotated_clusters, collection_dir, threads, timing_log)
        pa.create_nucleotide_alignment(annotated_clusters, collection_dir)
    if alignment == "nucleotide":
        pa.create_nuc_file_for_each_cluster(samples, gene_to_cluster_name, pan_ref_list, collection_dir)
        pa.run_mafft_nucleotide_alignment(annotated_clusters, collection_dir, threads, timing_log)

    pa.create_core_gene_alignment(annotated_clusters, gene_dictionary, samples, collection_dir)
