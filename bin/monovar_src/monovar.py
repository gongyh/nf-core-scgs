#!/usr/bin/env python
"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import os
from pickle import TRUE
import sys
import argparse
if sys.version_info.major == 3:
    import copyreg as copy_reg
else:
    import copy_reg
import types
import multiprocessing as mp
import numpy as np
from functools import partial
from datetime import datetime

import utils as U
import mp_genotype as M
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from calc_variant_prob import Calc_Var_Prob
from hzvcf import VCFDocument

# Required for using multiprocessing


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


copy_reg.pickle(types.MethodType, _pickle_method)

# Default values
Base_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}


def parse_args():
    parser = argparse.ArgumentParser(
        prog='MonoVar',
        usage='samtools mpileup --fasta-ref <REF> --bam-list <BAMS> [options] '
        '| monovar.py -b <BAMS> -f <REF> -o <OUTPUT> [options]',
        description='*** SNV calling on single-cell DNA data. ***'
    )
    parser.add_argument('--version', action='version', version='1.0.0_NB')

    parser.add_argument('-i', '--pileup', type=str, default='',
                        help='Pileup file. If not given, input is read from stdin.')
    parser.add_argument('-s', '--samples', type=str, default='',
                        help='File containing sample names, 1 sample per line. Required if no '
                        'bam files are provided.')
    parser.add_argument('-b', '--bam_file_list', type=str, default='',
                        help='List of Bam files in a text format.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output file (should end on ".vcf").')
    parser.add_argument('-f', '--ref_file', type=str, default='',
                        help='Reference genome file in .fa format.')
    parser.add_argument('-p', '--pe', type=float, default=0.002,
                        help='Probability of an error. Default = 0.002.')
    parser.add_argument('-a', '--pad', type=float, default=0.2,
                        help='Probability of an allelic dropout. Default = 0.2.')
    parser.add_argument('-t', '--threshold', type=float, default=0.05,
                        help='Threshold to use for calling variant. Default = 0.05.')
    parser.add_argument('-c', '--CF_flag', type=int, choices=[0, 1], default=1,
                        help='Flag for consensus filtering. Default = 1.')
    parser.add_argument('-m', '--cpus', type=int, default=1,
                        help='Number of cpus to use for threading. Default = 1.')
    # newly added arguments:
    parser.add_argument('-mpd', '--max_pileup_depth', type=int, default=10000,
                        help='Maximum pileup depth to take into account. Default = 10000.')
    parser.add_argument('-mrd', '--min_read_depth', type=int, default=1,
                        help='Minimum read depth required for SNV calling (per cell-locus). '
                        'Default = 1.')
    parser.add_argument('-th', '--theta', type=float, default=0.001,
                        help='Heterozygosity rate theta. Default = 0.001.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Turn of threading. Default = False.')

    args = parser.parse_args()
    return args


def main(args):
    if args.bam_file_list:
        # Obtain the RG IDs from the bam files
        with open(args.bam_file_list, 'r') as f:
            f_bam_list = f.read().strip().split('\n')
        bam_id_list = [U.get_BAM_RG(i.strip()) for i in f_bam_list]
    else:
        if not os.path.exists(args.samples):
            raise IOError('If input is mpileup file, a sample name file is '
                          'required! (-s|--sample <FILE>)')

        with open(args.samples, 'r') as f:
            bam_id_list = f.read().strip().split('\n')

    n_cells = len(bam_id_list)
    n_cells_threshold = n_cells / 2
    # no of possible alternate alleles {0, 1, 2, ..., 2m}
    max_allele_cnt = 2 * n_cells + 1
    # Table for all the required nCr
    nCr_matrix = U.get_nCr_mat(max_allele_cnt)
    # Dictionary for holding all the priors for different values of n
    prior_variant_dict = {i: U.calc_prior(args.theta, i, 1)
                          for i in range(n_cells + 1)}
    # Lists for all single_cell_ftr_pos objects, cells containing read support,
    # and cell containing alternate allele support
    all_single_cell_ftrs_list = np.zeros(n_cells, dtype=object)
    read_flag_row = np.zeros(n_cells, dtype=bool)
    alt_allele_flag_row = np.zeros(n_cells, dtype=bool)

    contigs = set([])
    # Initialize the pool of multiprocessing
    pool = mp.Pool(processes=args.cpus)

    # Open VCF file and print header
    vcf = VCFDocument(args.output, bam_id_list, args.ref_file)

    if args.pileup:
        with open(args.pileup, 'r') as f:
            lines = f.read().strip().split('\n')
        in_type = 'mpileup file'
    else:
        lines = sys.stdin
        in_type = 'bam stdin'

    if args.debug:
        print('\tStart iterating over {} with {} cells'
              .format(in_type, len(bam_id_list)))

    for line in lines:
        row = line.strip().split('\t')
        if line == '':
            continue

        contig = row[0]
        pos = int(row[1])
        refBase = row[2].strip().upper()
        original_refBase = refBase

        if refBase not in ['A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W']:
            continue

        total_depth = 0
        total_ref_depth = 0
        for i in range(1, n_cells + 1):
            if original_refBase in ['R', 'Y', 'M', 'K', 'S', 'W']:
                if int(row[3*i]) != 0:
                    total_bases = row[3*i + 1]
                    total_ins_del_rmvd_bases = U.ins_del_rmvd_original_bases(
                        total_bases)
                    total_start_and_end_bases = U.get_start_and_end(
                        total_ins_del_rmvd_bases)
                    total_count = U.get_base_count(total_start_and_end_bases)
                    total_count_descend_index = total_count.argsort()[::-1]
                    refBase = U.alt_deg_ref(
                        original_refBase, total_count_descend_index)
                else:
                    refBase = original_refBase

            curr_cell_pos_ftrs = Single_Cell_Ftrs_Pos(
                refBase, original_refBase, row[3*i: 3*i + 3])
            total_depth += curr_cell_pos_ftrs.depth
            total_ref_depth += curr_cell_pos_ftrs.refDepth
            all_single_cell_ftrs_list[i - 1] = curr_cell_pos_ftrs

        if total_depth <= 10 or total_depth == total_ref_depth:
            continue

        alt_count = total_depth - total_ref_depth
        alt_freq = float(alt_count) / total_depth

        # No reads supporting alternate allele, so no operations needed
        if alt_freq <= 0.01:
            continue
        # Cases that are to be prefiltered
        if total_depth > 30 and (alt_count <= 2 or alt_freq <= 0.001):
            continue

        # List for storing the sngl_cell_objs that have read support, will be
        # further used in the model
        read_supported_cell_list = []
        # Gloabal list for storing the alternate allele counts
        total_alt_allele_count = np.zeros(4, dtype=int)

        # Traverse through all the sngl_cell_ftr_obj and if has read support
        #   further calculate the other quantities
        for j, sngl_cell_ftr_obj in enumerate(all_single_cell_ftrs_list):
            read_flag = sngl_cell_ftr_obj.depth >= args.min_read_depth
            if read_flag:
                sngl_cell_ftr_obj.get_base_call_string_nd_quals()
                if original_refBase in ['R', 'Y', 'M', 'K', 'S', 'W']:
                    total_alt_allele_judgement = np.zeros(4, dtype=int)
                    total_alt_allele_judgement += sngl_cell_ftr_obj \
                        .get_deg_alt_allele_count(original_refBase)
                    if np.all(total_alt_allele_judgement == 0):
                        alt_allele_flag = False
                    else:
                        alt_allele_flag = True
                else:
                    alt_allele_flag = \
                        sngl_cell_ftr_obj.depth - sngl_cell_ftr_obj.refDepth != 0
                if alt_allele_flag:
                    # Update the list of total_alt_allele_count
                    if original_refBase in ['R', 'Y', 'M', 'K', 'S', 'W']:
                        total_alt_allele_count += sngl_cell_ftr_obj \
                            .get_deg_alt_allele_count(original_refBase)
                    else:
                        total_alt_allele_count += sngl_cell_ftr_obj \
                            .get_Alt_Allele_Count()
                # Populate the list of read supported cells
                read_supported_cell_list.append(sngl_cell_ftr_obj)
            else:
                alt_allele_flag = False
            read_flag_row[j] = read_flag
            alt_allele_flag_row[j] = alt_allele_flag

        # Operations on the single cells with read support
        # Number of cells with read support
        read_supported_n_cells = len(read_supported_cell_list)
        if read_supported_n_cells == 0:
            continue
        # Number of cells having read support
        read_smpl_count = read_flag_row.sum()
        # Number of cells having alternate allele support
        alt_smpl_count = alt_allele_flag_row.sum()
        # Update alt count
        alt_count = total_alt_allele_count.max()
        # Get the altBase
        if alt_count == 0:
            continue

        altBase = Base_dict[total_alt_allele_count.argmax()]

        # Calculate prior_allele_mat
        prior_allele_mat = U.get_prior_allele_mat(read_smpl_count, alt_smpl_count,
                                                  n_cells_threshold, total_depth, alt_freq, args.pe)

        # Get prior_variant_number distribution (Eq. 11)
        prior_var_no = prior_variant_dict[read_supported_n_cells]

        for cell in read_supported_cell_list:
            cell.store_addl_info(altBase, alt_freq, prior_allele_mat)

        # Obtain the value of probability of SNV
        var_prob_obj = Calc_Var_Prob(read_supported_cell_list)
        zero_var_prob, denominator = var_prob_obj \
            .calc_zero_var_prob(n_cells, args.max_pileup_depth, nCr_matrix,
                                args.pad, prior_var_no)

        # Skip of probability of SNV does not pass the threshold
        if zero_var_prob > args.threshold:
            continue

        # Global list for storing the indices of the cells having read support
        func = partial(M.get_info_string, read_supported_cell_list, n_cells,
                       nCr_matrix, prior_var_no, denominator)

        if args.debug:
            output = [func(i) for i in range(read_supported_n_cells)]
        else:
            output = pool.map(func, range(read_supported_n_cells))

        barcode = []
        info_list = []
        for single_cell_ftrs_list in all_single_cell_ftrs_list:
            if single_cell_ftrs_list.depth == 0:
                info_list.append('./.')
                barcode.append('X')
            else:
                info_cell, barcode_cell = output.pop(0)
                info_list.append(info_cell)
                barcode.append(barcode_cell)

        if zero_var_prob == 0:
            qual = 99
        else:
            qual = min(99, int(np.round(-10 * np.log10(zero_var_prob))))

        AC, AF, AN = U.calc_chr_count(barcode)
        if total_ref_depth > 0:
            baseQranksum = U.calc_base_q_rank_sum(read_supported_cell_list)
        else:
            baseQranksum = 0
        QD = U.calc_qual_depth(barcode, all_single_cell_ftrs_list, qual)
        SOR = U.calc_strand_bias(read_supported_cell_list, alt_count)
        max_prob_ratio = U.find_max_prob_ratio(var_prob_obj.matrix)
        PSARR = U.calc_per_smpl_alt_ref_ratio(
            total_ref_depth, alt_count, read_smpl_count, alt_smpl_count)

        # Write record/line to vcf output
        info_str = 'AC={ac};AF={af:.2f};AN={an};BaseQRankSum={bqrs:.2f};DP={dp};' \
            'QD={qd:.4f};SOR={sor:.2f};MPR={mpr:.2f};PSARR={psarr:.2f}' \
            .format(ac=AC, af=AF, an=AN, bqrs=baseQranksum, dp=total_depth,
                    qd=QD, sor=SOR, mpr=max_prob_ratio, psarr=PSARR)
        sample_str = '\t'.join(info_list)

        if args.CF_flag:
            if U.consensus_filter(barcode):
                filter_str = 'PASS'
            else:
                filter_str = 'NoConsensus'
        else:
            filter_str = '.'

        vcf_rec_data = [contig, str(pos), '.', original_refBase, altBase, str(qual),
                        filter_str, info_str, 'GT:AD:DP:GQ:PL', sample_str]

        vcf.append_record(vcf_rec_data)
        contigs.add(contig)

    if args.debug:
        print('\tCreating final vcf file {}'.format(args.output))

    vcf.close_records()
    vcf.add_contigs(contigs)
    vcf.add_header()


if __name__ == '__main__':
    args = parse_args()
    print('Start Monovar_NB: {:%Y%m%d_%H:%M:%S}'.format(datetime.now()))
    main(args)
    print('Stop  Monovar_NB: {:%Y%m%d_%H:%M:%S}'.format(datetime.now()))
