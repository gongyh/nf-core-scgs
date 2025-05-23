# -*- coding: utf-8 -*-
"""
Wrapper for pipelines
"""
import os
import subprocess
import logging

from pan_genome import data_preparation
from pan_genome import main_pipeline
from pan_genome import add_sample_pipeline
from pan_genome import annotate
from pan_genome import alignment
from pan_genome import output
from pan_genome import post_analysis

logger = logging.getLogger(__name__)


def run_init_pipeline(samples, collection_dir, temp_dir, baseDir, args, timing_log, resume):
    """
    Run initial pan-genome analysis.

    Pipeline description:
        (1) extract gene sequences
        (2) cluster by CD-HIT
        (3) All-agaist-all comparision by BLASTP
        (4) cluster by Markov clustering algorithms - MCL
        (5) annotate gene clusters
        (6) Create multiple sequences alignment by abPOA
        (7) write output

    Parameters
    ----------
    samples : list
        list of samples
    collection_dir : path
        collection directory
    temp_dir : path
        temporary directory
    baseDir : path
        directory of panta
    args : object
        Command-line input arguments
    timing_log : path
        path of time.log
    resume : list
        A boolean inside a list
        If True, resume previous analysis
    """
    gene_dictionary, gene_position = data_preparation.extract_proteins(samples, collection_dir, args)

    combined_faa = data_preparation.combine_proteins(
        collection_dir=collection_dir, out_dir=temp_dir, samples=samples, timing_log=timing_log, resume=resume
    )

    cd_hit_represent_fasta, cd_hit_clusters = main_pipeline.run_cd_hit(
        faa_file=combined_faa, out_dir=temp_dir, threads=args.threads, timing_log=timing_log, resume=resume
    )

    blast_result = main_pipeline.pairwise_alignment(
        diamond=args.diamond,
        database_fasta=cd_hit_represent_fasta,
        query_fasta=cd_hit_represent_fasta,
        out_dir=os.path.join(temp_dir, "blast"),
        timing_log=timing_log,
        evalue=args.evalue,
        max_target_seqs=2000,
        identity=args.identity,
        query_coverage=args.coverage,
        threads=args.threads,
        resume=resume,
    )

    mcl_file = main_pipeline.cluster_with_mcl(
        out_dir=temp_dir, blast_result=blast_result, timing_log=timing_log, resume=resume
    )

    unsplit_clusters = main_pipeline.reinflate_clusters(cd_hit_clusters=cd_hit_clusters, mcl_file=mcl_file)

    clusters = post_analysis.split_paralogs(gene_dictionary, gene_position, unsplit_clusters, args.split)

    ## TODO: call variants from msa, output vcf file
    if args.fasta == None:  # if input data is annotated genomes
        clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=clusters, gene_dictionary=gene_dictionary
        )
        output.output_cluster_info(clusters, clusters_annotation, gene_dictionary, samples, collection_dir)
        representative_fasta = alignment.create_msa_init_pipeline(
            clusters, samples, collection_dir, baseDir, args.threads
        )
    else:  # if input data is genome assembly
        representative_fasta = alignment.create_msa_init_pipeline(
            clusters, samples, collection_dir, baseDir, args.threads
        )
        # clusters_annotation = annotate.annotate_cluster_fasta(
        #    unlabeled_clusters=clusters,
        #    rep_fasta = representative_fasta,
        #    temp_dir=temp_dir,
        #    baseDir = baseDir,
        #    timing_log=timing_log,
        #    threads = args.threads)
        output.output_cluster_info(clusters, [], gene_dictionary, samples, collection_dir)

    # create a new gene info file
    gene_file = os.path.join(collection_dir, "gene_info.tsv")
    if os.path.exists(gene_file):
        os.remove(gene_file)  # remove existing file to create a new one
    output.output_gene_info(clusters, gene_dictionary, collection_dir)


def run_add_pipeline(
    new_samples, old_represent_faa, previous_clusters, collection_dir, temp_dir, baseDir, args, timing_log, resume
):
    """
    Add new samples to previous collection.

    Pipeline description:
        (1) extract gene sequences of new samples
        (2) match new sequence with previous clusters by CD-HIT-2D
        (3) cluster the remain sequences by CD-HIT
        (4) match new sequence with previous clusters by BLASTP
        (5) All-agaist-all comparision by BLASTP
        (6) cluster by Markov clustering algorithms - MCL
        (7) annotate gene clusters
        (8) Create multiple sequences alignment by abPOA
        (9) update output

    Parameters
    ----------
    new_samples : list
        list of new samples
    old_represent_faa : path
        path of reference_pangenome.fasta
    previous_clusters : list  of []
        each empty list correspond to a previous cluster
    collection_dir : path
        collection directory
    temp_dir : path
        temporary directory
    baseDir : path
        directory of panta
    args : object
        Command-line input arguments
    timing_log : path
        path of time.log
    resume : list
        A boolean inside a list
        If True, resume previous analysis
    """
    gene_dictionary, _ = data_preparation.extract_proteins(new_samples, collection_dir, args)

    new_combined_faa = data_preparation.combine_proteins(
        collection_dir=collection_dir, out_dir=temp_dir, samples=new_samples, timing_log=timing_log, resume=resume
    )

    notmatched_faa, cd_hit_2d_clusters = add_sample_pipeline.run_cd_hit_2d(
        database_1=old_represent_faa,
        database_2=new_combined_faa,
        out_dir=temp_dir,
        threads=args.threads,
        timing_log=timing_log,
        resume=resume,
    )

    add_sample_pipeline.add_gene_cd_hit_2d(previous_clusters, cd_hit_2d_clusters)

    num_seq = subprocess.run(f'grep ">" {notmatched_faa} | wc -l', capture_output=True, text=True, shell=True)
    if int(num_seq.stdout.rstrip()) == 0:
        new_clusters = []
    else:
        notmatched_represent_faa, cd_hit_clusters = main_pipeline.run_cd_hit(
            faa_file=notmatched_faa, out_dir=temp_dir, threads=args.threads, timing_log=timing_log, resume=resume
        )

        blast_1_result = main_pipeline.pairwise_alignment(
            diamond=args.diamond,
            database_fasta=old_represent_faa,
            query_fasta=notmatched_represent_faa,
            out_dir=os.path.join(temp_dir, "blast1"),
            timing_log=timing_log,
            evalue=args.evalue,
            max_target_seqs=2000,
            identity=args.identity,
            query_coverage=args.coverage,
            threads=args.threads,
            resume=resume,
        )

        remain_fasta = add_sample_pipeline.add_gene_blast(
            previous_clusters=previous_clusters,
            cd_hit_clusters=cd_hit_clusters,
            blast_result=blast_1_result,
            fasta_file=notmatched_represent_faa,
            out_dir=temp_dir,
        )

        num_seq = subprocess.run(f'grep ">" {remain_fasta} | wc -l', capture_output=True, text=True, shell=True)
        if int(num_seq.stdout.rstrip()) == 0:
            new_clusters = []
        else:
            blast_2_result = main_pipeline.pairwise_alignment(
                diamond=args.diamond,
                database_fasta=remain_fasta,
                query_fasta=remain_fasta,
                out_dir=os.path.join(temp_dir, "blast2"),
                timing_log=timing_log,
                evalue=args.evalue,
                max_target_seqs=2000,
                identity=args.identity,
                query_coverage=args.coverage,
                threads=args.threads,
                resume=resume,
            )

            mcl_file = main_pipeline.cluster_with_mcl(
                out_dir=temp_dir, blast_result=blast_2_result, timing_log=timing_log, resume=resume
            )

            new_clusters = main_pipeline.reinflate_clusters(cd_hit_clusters=cd_hit_clusters, mcl_file=mcl_file)

    # annotate clusters, create gene alignment and output
    ## TODO: reannotate old clusters, which could be changed after
    # adding new sequences
    if args.fasta == None:  # if input data is annotated genome
        # annotate new clusters only
        new_clusters_annotation = annotate.annotate_cluster_gff(
            unlabeled_clusters=new_clusters, gene_dictionary=gene_dictionary
        )
        # update result files
        output.update_cluster_info(
            previous_clusters,
            new_clusters,
            new_clusters_annotation,
            gene_dictionary,
            new_samples,
            temp_dir,
            collection_dir,
        )
        # add new seq to previous msa, create msa for new clusters
        new_representative_fasta = alignment.create_msa_add_pipeline(
            previous_clusters, new_clusters, new_samples, collection_dir, baseDir, args.threads
        )
    else:  # if input data is genome assembly
        # add new seq to previous msa, create msa for new clusters
        new_representative_fasta = alignment.create_msa_add_pipeline(
            previous_clusters, new_clusters, new_samples, collection_dir, baseDir, args.threads
        )
        # annotate new clusters only
        new_clusters_annotation = annotate.annotate_cluster_fasta(
            unlabeled_clusters=new_clusters,
            rep_fasta=new_representative_fasta,
            temp_dir=temp_dir,
            baseDir=baseDir,
            timing_log=timing_log,
            threads=args.threads,
        )
        # update result files
        output.update_cluster_info(
            previous_clusters,
            new_clusters,
            new_clusters_annotation,
            gene_dictionary,
            new_samples,
            temp_dir,
            collection_dir,
        )

    # append new gene to previous gene info file
    previous_clusters.extend(new_clusters)
    output.output_gene_info(previous_clusters, gene_dictionary, collection_dir)
