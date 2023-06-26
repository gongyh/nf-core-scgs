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

from asyncore import file_dispatcher
import copy
from curses.ascii import alt
import re
import numpy as np
import pysam
from scipy import stats
from base_q_ascii import base_q_dict, base_q_int_dict
from alleles_prior import get_prior_matrix
import random

Base_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}


def get_nCr_mat(max_allele_cnt):
    factorial_list = [np.math.factorial(i) for i in range(max_allele_cnt)]

    ncr_mat = np.zeros((max_allele_cnt, max_allele_cnt))
    for i in range(max_allele_cnt):
        for j in range(max_allele_cnt):
            ncr_mat[j, i] = factorial_list[j] / \
                (factorial_list[i] * factorial_list[j - i])
    return ncr_mat


def get_ref_count(read_base, ref):
    forward_ref_c = read_base.count('.')
    reverse_ref_c = read_base.count(',')
    ref_count = forward_ref_c + reverse_ref_c
    return forward_ref_c, reverse_ref_c, ref_count


def get_deg_ref_count(read_base, ref):
    if ref == 'R':
        forward_ref_c = read_base.count('A') + read_base.count('G')
        reverse_ref_c = read_base.count('a') + read_base.count('g')
    elif ref == 'Y':
        forward_ref_c = read_base.count('C') + read_base.count('T')
        reverse_ref_c = read_base.count('c') + read_base.count('t')
    elif ref == 'M':
        forward_ref_c = read_base.count('A') + read_base.count('C')
        reverse_ref_c = read_base.count('a') + read_base.count('c')
    elif ref == 'K':
        forward_ref_c = read_base.count('G') + read_base.count('T')
        reverse_ref_c = read_base.count('g') + read_base.count('t')
    elif ref == 'S':
        forward_ref_c = read_base.count('G') + read_base.count('C')
        reverse_ref_c = read_base.count('g') + read_base.count('c')
    else:
        forward_ref_c = read_base.count('A') + read_base.count('T')
        reverse_ref_c = read_base.count('a') + read_base.count('t')
    ref_count = forward_ref_c + reverse_ref_c
    return forward_ref_c, reverse_ref_c, ref_count


def copy_list_but_one(i_list, index):
    nu_list = [i_list[i] for i in range(index)]
    for i in range(index + 1, len(i_list)):
        nu_list.append(i_list[i])
    return nu_list


def find_indel(string, pattern):
    """ Finds all the occurrences of pattern in string. Removes all 
    occurrences of pattern from the string and returns list of found patterns
    and the new string with patterns removed
    """
    l = [x.group() for x in re.finditer(pattern, string)]
    len_indel = [int(re.split(r'(\d+)', i)[1])
                 for i in l]  # Get the lengths of the patterns
    spans = [i.span() for i in re.finditer(pattern, string)]  # Get the spans
    newspan = []  # Change the spans according to the integer length of the pattern
    for i in range(len(len_indel)):
        new_end = spans[i][0] + 1 + len_indel[i] + len(str(len_indel[i]))
        newspan.append((spans[i][0], new_end))
    # Find the instances of the patterns
    final_indel_list = [string[i1:i2] for (i1, i2) in newspan]
    new_string = string
    for i in final_indel_list:
        new_string = new_string.replace(i, '')
    return final_indel_list, new_string


def get_alt_count(string):
    ins_list, ins_rmvd_str = find_indel(string, '\+[0-9]+[ACGTNacgtn]+')
    del_list, del_ins_rmvd_str = \
        find_indel(ins_rmvd_str, '-[0-9]+[ACGTNacgtn]+')
    ins_count = len(ins_list)
    del_count = len(del_list)
    A_cnt = del_ins_rmvd_str.count('A') + del_ins_rmvd_str.count('a')
    T_cnt = del_ins_rmvd_str.count('T') + del_ins_rmvd_str.count('t')
    G_cnt = del_ins_rmvd_str.count('G') + del_ins_rmvd_str.count('g')
    C_cnt = del_ins_rmvd_str.count('C') + del_ins_rmvd_str.count('c')
    N_cnt = del_ins_rmvd_str.count('N') + del_ins_rmvd_str.count('n')
    return del_ins_rmvd_str, ins_count, del_count, A_cnt, T_cnt, G_cnt, \
        C_cnt, N_cnt


def get_count_start_and_end(s):
    end_counts = s.count('$')
    ns = s.replace('$', '')
    start_counts = 0
    i = 0
    fs = ''
    while (i < len(ns)):
        if ns[i] == '^':
            i += 2
            start_counts += 1
        else:
            fs = fs + ns[i]
            i += 1
    return start_counts, end_counts, fs


def get_base_call_string(s, ref):
    """ Removes unwanted characters from s and then replace . and , with ref, 
        finally returns a string that contains all the observed bases
    """
    l = {'.': ref, ',': ref, '*': ref, 'a': 'A', 'A': 'A', 'c': 'C', 'C': 'C',
         't': 'T', 'T': 'T', 'g': 'G', 'G': 'G'}
    sn = ''
    for i in s:
        try:
            sn += l[i]
        except KeyError:
            pass
    return sn


def get_base_qual_list(qual_str):
    """ Returns the base quality scores as a list """
    len_s = len(qual_str)
    base_q_list = np.zeros(len_s)
    base_q_int_list = np.zeros(len_s)
    for i, qual_base in enumerate(qual_str):
        base_q_list[i] = base_q_dict[qual_base]
        base_q_int_list[i] = base_q_int_dict[qual_base]
    return base_q_list, base_q_int_list


def calc_strand_bias(cell_ftr_pos_list, Alt_count):
    f_ref_count = 0
    f_alt_count = 0
    r_ref_count = 0
    r_alt_count = 0
    for cell_ftr_pos in cell_ftr_pos_list:
        fr, fa, rr, ra = cell_ftr_pos.get_strand_bias_info()
        f_ref_count += fr
        f_alt_count += fa
        r_ref_count += rr
        r_alt_count += ra

    if f_ref_count == 0 and r_ref_count == 0:
        oddsRatio = 0.0
    else:
        cont_table = np.array(
            [[f_ref_count, r_ref_count], [f_alt_count, r_alt_count]]
        )
        oddsRatio, pval = stats.fisher_exact(cont_table)

    return oddsRatio


def calc_prior(theta, n_cells, flag):
    prior_variant_number = []
    if flag == 1:
        for i in range(0, 2 * n_cells + 1):
            if ((i == 0) | (i == 2 * n_cells)):
                lst = [1.0 / i for i in range(1, 2 * n_cells)]
                prob = 0.5 * (1 - (theta * sum(lst)))
            else:
                prob = theta / i
            prior_variant_number.append(prob)
    elif flag == 2:
        norm_const = 0
        for i in range(1, 2 * n_cells):
            if (i == n_cells):
                norm_const += 2 * n_cells
            else:
                norm_const += float(n_cells) / abs(n_cells - i)
        for i in range(1, 2 * n_cells):
            if (i == n_cells):
                prob = 2 * n_cells * theta / norm_const
            else:
                prob = ((n_cells * theta) / abs(n_cells - i)) / norm_const
            prior_variant_number.append(prob)
        sp = sum(prior_variant_number)
        p_l0 = 0.5 * (1 - sp)
        p_l2n = 0.5 * (1 - sp)
        prior_variant_number.insert(0, p_l0)
        prior_variant_number.append(p_l2n)
    elif flag == 3:
        for i in range(0, 2 * n_cells + 1):
            prob = 1. / 2 * n_cells + 1
            prior_variant_number.append(prob)
    return prior_variant_number


def find_max_prob_ratio(matrix):
    l_0_prob = matrix.denom_prob_matrix[0, -1]
    if l_0_prob == 0:
        subtracting_max_prob = -743.7469
    else:
        subtracting_max_prob = np.log(l_0_prob)

    allele_count = matrix.denom_prob_matrix[:, -1].argmax()
    max_prob = matrix.denom_prob_matrix[allele_count, -1]

    if max_prob <= 0:
        log_max_prob = -743.7469
    else:
        log_max_prob = np.log(max_prob)
    max_prob_ratio = log_max_prob - subtracting_max_prob
    return max_prob_ratio


def get_prior_allele_mat(read_smpl_count, alt_smpl_count,
                         cell_no_threshold, total_depth, Alt_freq, pe):
    if read_smpl_count > cell_no_threshold - 1 and alt_smpl_count == 1:
        prior_mat = get_prior_matrix(0.2)
    elif read_smpl_count > cell_no_threshold and alt_smpl_count == 2 \
            and total_depth > 30 and Alt_freq < 0.1:
        prior_mat = get_prior_matrix(0.1)
    else:
        prior_mat = get_prior_matrix(pe)
    return prior_mat


def calc_chr_count(barcode):
    AC = 0
    AN = 0
    for c in barcode:
        if c == 'X':
            continue
        else:
            AN += 2
            AC += c
    AF = float(AC) / AN
    return AC, AF, AN


def calc_base_q_rank_sum(read_supported_cell_list):
    ref_list = []
    alt_list = []
    for cell_ftr_info in read_supported_cell_list:
        for i, base in enumerate(cell_ftr_info.final_bases):
            if base == cell_ftr_info.refBase:
                ref_list.append(cell_ftr_info.base_qual_int_val_list[i])
            elif base == cell_ftr_info.altBase:
                alt_list.append(cell_ftr_info.base_qual_int_val_list[i])

    baseQranksum, pVal = stats.ranksums(alt_list, ref_list)
    return baseQranksum


def calc_qual_depth(barcode, all_single_cell_ftrs_list, qual):
    depth = 0
    for i, c in enumerate(barcode):
        if c == 'X' or c == 0:
            continue
        depth += all_single_cell_ftrs_list[i].depth
    if depth > 0:
        qual_depth = float(qual) / depth
    else:
        qual_depth = qual
    return qual_depth


def get_BAM_RG(bam_file):
    rows = pysam.view("-H", bam_file)
    for r in rows:
        if r.startswith('@RG'):
            r_l = r.split('\t')
            id_list = r_l[1].split(':')
            return id_list[1]

    bam_id = bam_file.split('/')[-1].strip('..').strip('~')
    return bam_id


def calc_per_smpl_alt_ref_ratio(total_ref_depth, alt_count,
                                read_smpl_count, alt_smpl_count):
    if total_ref_depth == 0:
        denom = 1
    else:
        denom = float(total_ref_depth) / read_smpl_count
    return (float(alt_count) / alt_smpl_count) / denom


def consensus_filter(barcode):
    g_count = 0
    for c in barcode:
        if c == 1 or c == 2:
            g_count += 1

    if g_count > 1:
        return True
    else:
        return False


def ins_del_rmvd_original_bases(original_bases):
    if original_bases.count('+') + original_bases.count('-') == 0:
        ins_del_rmvd_bases = original_bases
    else:
        ins_list = []
        del_list = []
        cp_original_bases = original_bases
        ins_list, ins_rmvd_bases = \
            find_indel(cp_original_bases,
                       '\+[0-9]+[ACGTNRYMKSWBVDHZacgtnrymkswbvdhz]+')
        del_list, ins_del_rmvd_bases = \
            find_indel(ins_rmvd_bases,
                       '-[0-9]+[ACGTNRYMKSWBVDHZacgtnrymkswbvdhz]+')
    return ins_del_rmvd_bases


def get_start_and_end(s):
    ns = s.replace('$', '')
    i = 0
    fs = ''
    while (i < len(ns)):
        if ns[i] == '^':
            i += 2
        else:
            fs = fs + ns[i]
            i += 1
    return fs


def get_base_count(final_bases):
    A_cnt = final_bases.count('A') \
        + final_bases.count('a')
    T_cnt = final_bases.count('T') \
        + final_bases.count('t')
    G_cnt = final_bases.count('G') \
        + final_bases.count('g')
    C_cnt = final_bases.count('C') \
        + final_bases.count('c')
    return np.array([A_cnt, T_cnt, G_cnt, C_cnt], dtype=int)


def alt_deg_ref(refBase, total_count_descend_index):
    if refBase == 'R':
        for i in range(4):
            if total_count_descend_index[i] == 0 or total_count_descend_index[i] == 2:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    elif refBase == 'Y':
        for i in range(4):
            if total_count_descend_index[i] == 1 or total_count_descend_index[i] == 3:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    elif refBase == 'M':
        for i in range(4):
            if total_count_descend_index[i] == 0 or total_count_descend_index[i] == 3:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    elif refBase == 'K':
        for i in range(4):
            if total_count_descend_index[i] == 1 or total_count_descend_index[i] == 2:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    elif refBase == 'S':
        for i in range(4):
            if total_count_descend_index[i] == 2 or total_count_descend_index[i] == 3:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    else:
        for i in range(4):
            if total_count_descend_index[i] == 0 or total_count_descend_index[i] == 1:
                ref_alt_index = total_count_descend_index[i]
                break
        ref_alt = Base_dict[ref_alt_index]
    return ref_alt
