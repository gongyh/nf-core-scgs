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

import numpy as np
import utils as U

G_MAP = {'CA': 'AC', 'GA': 'AG', 'TA': 'AT',
         'GC': 'CG', 'TC': 'CT', 'TG': 'GT'}


class Single_Cell_Ftrs_Pos:

    # Constructor takes the list of info for the current position as input
    # ([contig, loc, ref, depth, primary_bases, base_q])
    def __init__(self, refBase, ori_refBase, current_pos_info_list):
        if (int(current_pos_info_list[0]) == 0):
            self.depth = 0
            self.refDepth = 0
        else:
            self.refBase = refBase
            self.ori_refBase = ori_refBase
            self.depth = int(current_pos_info_list[0])
            self.primary_bases = current_pos_info_list[1]
            self.base_q = current_pos_info_list[2]

            if ori_refBase in ['R', 'Y', 'M', 'K', 'S', 'W']:
                self.forward_ref_count, self.reverse_ref_count, self.refDepth = \
                    U.get_deg_ref_count(self.primary_bases, self.ori_refBase)
            else:
                self.forward_ref_count, self.reverse_ref_count, self.refDepth = \
                    U.get_ref_count(self.primary_bases, self.refBase)

    def __str__(self):
        try:
            out_str = 'refBase={},depth={},refDepth={}\n' \
                'primary_bases=\t{}\nbase_q=\t\t{}\n' \
                .format(self.refBase, self.depth, self.refDepth,
                        self.primary_bases, self.base_q)
        except AttributeError:
            out_str = 'depth={},refDepth={}'.format(self.depth, self.refDepth)
        return out_str

    # Remove the insertions and deletions from the primary_bases and also
    # count the number of insertions and deletions

    def get_ins_del_rmvd_bases(self):
        if self.primary_bases.count('+') + self.primary_bases.count('-') == 0:
            self.ins_count = 0
            self.del_count = 0
            self.ins_list = []
            self.del_list = []
            self.ins_del_rmvd_bases = self.primary_bases
        else:
            cp_primary_bases = self.primary_bases
            self.ins_list, ins_rmvd_bases = \
                U.find_indel(cp_primary_bases,
                             '\+[0-9]+[ACGTNRYMKSWBVDHZacgtnrymkswbvdhz]+')
            self.ins_count = len(self.ins_list)
            self.del_list, self.ins_del_rmvd_bases = \
                U.find_indel(ins_rmvd_bases,
                             '-[0-9]+[ACGTNRYMKSWBVDHZacgtnrymkswbvdhz]+')
            self.del_count = len(self.del_list)

    # Function that calculates the base calling errors from base qual scores

    def get_base_qual_vals(self):
        self.base_qual_val_list, self.base_qual_int_val_list = \
            U.get_base_qual_list(self.base_q)

    # After removal of insertions and deletions we create the base call string
    # to be used by the model

    def get_base_calls(self, ref):
        self.start_read_counts, self.end_read_counts, self.start_end_ins_del_rmvd_bases = \
            U.get_count_start_and_end(self.ins_del_rmvd_bases)
        self.final_bases = \
            U.get_base_call_string(self.start_end_ins_del_rmvd_bases, ref)

    def get_base_call_string_nd_quals(self):
        self.get_ins_del_rmvd_bases()
        self.get_base_qual_vals()
        self.get_base_calls(self.refBase)

    # Function that calculates the numbers of alternate alleles in the
    # ins_del_rmvd_bases

    def get_Alt_Allele_Count(self):
        A_cnt = self.start_end_ins_del_rmvd_bases.count('A') \
            + self.start_end_ins_del_rmvd_bases.count('a')
        T_cnt = self.start_end_ins_del_rmvd_bases.count('T') \
            + self.start_end_ins_del_rmvd_bases.count('t')
        G_cnt = self.start_end_ins_del_rmvd_bases.count('G') \
            + self.start_end_ins_del_rmvd_bases.count('g')
        C_cnt = self.start_end_ins_del_rmvd_bases.count('C') \
            + self.start_end_ins_del_rmvd_bases.count('c')
        return np.array([A_cnt, T_cnt, G_cnt, C_cnt], dtype=int)

    def get_deg_alt_allele_count(self, original_refBase):
        if original_refBase == 'R':
            A_cnt = 0
            T_cnt = self.start_end_ins_del_rmvd_bases.count('T') \
                + self.start_end_ins_del_rmvd_bases.count('t')
            G_cnt = 0
            C_cnt = self.start_end_ins_del_rmvd_bases.count('C') \
                + self.start_end_ins_del_rmvd_bases.count('c')
        elif original_refBase == 'Y':
            A_cnt = self.start_end_ins_del_rmvd_bases.count('A') \
                + self.start_end_ins_del_rmvd_bases.count('a')
            T_cnt = 0
            G_cnt = self.start_end_ins_del_rmvd_bases.count('G') \
                + self.start_end_ins_del_rmvd_bases.count('g')
            C_cnt = 0
        elif original_refBase == 'M':
            A_cnt = 0
            T_cnt = self.start_end_ins_del_rmvd_bases.count('T') \
                + self.start_end_ins_del_rmvd_bases.count('t')
            G_cnt = self.start_end_ins_del_rmvd_bases.count('G') \
                + self.start_end_ins_del_rmvd_bases.count('g')
            C_cnt = 0
        elif original_refBase == 'K':
            A_cnt = self.start_end_ins_del_rmvd_bases.count('A') \
                + self.start_end_ins_del_rmvd_bases.count('a')
            T_cnt = 0
            G_cnt = 0
            C_cnt = self.start_end_ins_del_rmvd_bases.count('C') \
                + self.start_end_ins_del_rmvd_bases.count('c')
        elif original_refBase == 'S':
            A_cnt = self.start_end_ins_del_rmvd_bases.count('A') \
                + self.start_end_ins_del_rmvd_bases.count('a')
            T_cnt = self.start_end_ins_del_rmvd_bases.count('T') \
                + self.start_end_ins_del_rmvd_bases.count('t')
            G_cnt = 0
            C_cnt = 0
        else:
            A_cnt = 0
            T_cnt = 0
            G_cnt = self.start_end_ins_del_rmvd_bases.count('G') \
                + self.start_end_ins_del_rmvd_bases.count('g')
            C_cnt = self.start_end_ins_del_rmvd_bases.count('C') \
                + self.start_end_ins_del_rmvd_bases.count('c')
        return np.array([A_cnt, T_cnt, G_cnt, C_cnt], dtype=int)

    # Store altBase, Alt_freq, prior_allele_mat
    def store_addl_info(self, altBase, Alt_freq, prior_allele_mat):
        self.altBase = altBase
        self.Alt_freq = Alt_freq
        self.forward_alt_count = self.start_end_ins_del_rmvd_bases \
            .count(self.altBase)
        self.reverse_alt_count = self.start_end_ins_del_rmvd_bases \
            .count(self.altBase.lower())
        self.alt_count = self.forward_alt_count + self.reverse_alt_count

        self.prior_allele_mat = prior_allele_mat

    def get_strand_bias_info(self):
        return (self.forward_ref_count, self.forward_alt_count,
                self.reverse_ref_count, self.reverse_alt_count)

    # Calculate probability of data given genotype gt

    def calc_prob_gt(self, gt, max_depth, start=0):
        if self.ori_refBase not in ['R', 'Y', 'M', 'K', 'S', 'W']:
            ub = min(len(self.base_qual_val_list),
                     len(self.final_bases), max_depth)
            val = np.ones(ub - start)
            for i in range(start, ub, 1):
                curr_base = self.final_bases[i]
                curr_base_key = (gt, curr_base)
                curr_err = self.base_qual_val_list[i]
                if curr_err > 0.5:
                    continue
                prob_i = self.prior_allele_mat[curr_base_key]
                # Eq. 1, Eq. 2, Eq. 4
                prob = curr_err * (1 - prob_i) / 3 + (1 - curr_err) * prob_i
                val[i - start] = prob
        else:
            if self.ori_refBase == 'R':
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'A' and gt[0] == 'G':
                        g = gt.replace('G', 'A')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'G' and gt[0] == 'A':
                        g = gt.replace('A', 'G')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
            elif self.ori_refBase == 'Y':
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'C' and gt[0] == 'T':
                        g = gt.replace('T', 'C')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'T' and gt[0] == 'C':
                        g = gt.replace('C', 'T')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
            elif self.ori_refBase == 'M':
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'A' and gt[0] == 'C':
                        g = gt.replace('C', 'A')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'C' and gt[0] == 'A':
                        g = gt.replace('A', 'C')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
            elif self.ori_refBase == 'K':
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'G' and gt[0] == 'T':
                        g = gt.replace('T', 'G')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'T' and gt[0] == 'G':
                        g = gt.replace('G', 'T')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
            elif self.ori_refBase == 'S':
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'G' and gt[0] == 'C':
                        g = gt.replace('C', 'G')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'C' and gt[0] == 'G':
                        g = gt.replace('G', 'C')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
            else:
                ub = min(len(self.base_qual_val_list),
                         len(self.final_bases), max_depth)
                val = np.ones(ub - start)
                for i in range(start, ub, 1):
                    curr_base = self.final_bases[i]
                    if curr_base == 'A' and gt[0] == 'T':
                        g = gt.replace('T', 'A')
                        curr_base_key = (g, curr_base)
                    if curr_base == 'T' and gt[0] == 'A':
                        g = gt.replace('A', 'T')
                        curr_base_key = (g, curr_base)
                    curr_base_key = (gt, curr_base)
                    curr_err = self.base_qual_val_list[i]
                    if curr_err > 0.5:
                        continue
                    prob_i = self.prior_allele_mat[curr_base_key]
                    # Eq. 1, Eq. 2, Eq. 4
                    prob = curr_err * (1 - prob_i) / 3 + \
                        (1 - curr_err) * prob_i
                    val[i - start] = prob
        return np.prod(val)

    # Function to calculate likelihood of hetero genotype for ub = all

    def Prob_Reads_Given_Genotype_hetero(self, g, ub, pad):
        prob_0_50d = self.cell_prob_0_50d
        prob_2_50d = self.cell_prob_2_50d
        prob_1_50d = self.calc_prob_gt(g, min(ub, 100))

        prob_50d = (1 - pad) * prob_1_50d \
            + pad / 2 * (prob_0_50d + prob_2_50d)
        self.cell_prob_1_50d = prob_50d
        if ub <= 100:
            return prob_50d
        else:
            prob_0 = self.cell_prob_0
            prob_2 = self.cell_prob_2
            # p(d|g=1, ADO = False)
            prob_1 = self.cell_prob_1_50d * self.calc_prob_gt(g, ub, 100)
            # Eq. 3: (1 - p_ad) * p(d|g=1, ADO=False) + p_ad * p(d|g=1, ADO=True)
            prob = (1 - pad) * prob_1 + pad / 2 * (prob_0 + prob_2)
            return prob

    def Prob_Reads_Given_Genotype(self, genotype_flag, max_depth, pad):
        if (self.altBase == ''):
            if (genotype_flag != 0):
                self.cell_prob_2 = 0
                self.cell_prob_1 = 0
                return 0.0
            else:
                g = self.refBase + self.refBase
        else:
            if (genotype_flag == 0):
                g = self.refBase + self.refBase
            elif (genotype_flag == 2):
                g = self.altBase + self.altBase
            else:
                g = self.refBase + self.altBase
                g = G_MAP.get(g, g)

        ub = min(len(self.base_qual_val_list),
                 len(self.final_bases), max_depth)

        if (genotype_flag == 0):
            self.cell_prob_0_50d = self.calc_prob_gt(g, min(100, ub))
            if ub <= 100:
                self.cell_prob_0 = self.cell_prob_0_50d
            else:
                self.cell_prob_0 = self.cell_prob_0_50d \
                    * self.calc_prob_gt(g, ub, 100)
            prob = self.cell_prob_0
        elif (genotype_flag == 2):
            self.cell_prob_2_50d = self.calc_prob_gt(g, min(100, ub))
            if ub <= 100:
                self.cell_prob_2 = self.cell_prob_2_50d
            else:
                self.cell_prob_2 = self.cell_prob_2_50d \
                    * self.calc_prob_gt(g, ub, 100)
            prob = self.cell_prob_2
        elif (genotype_flag == 1):
            self.cell_prob_1 = self.Prob_Reads_Given_Genotype_hetero(
                g, ub, pad)
            prob = self.cell_prob_1

        return prob

    def Prob_Reads_Given_Genotype_50d(self, gt_flag):
        if gt_flag == 0:
            self.cell_prob_0 = self.cell_prob_0_50d
            return self.cell_prob_0_50d
        elif gt_flag == 1:
            self.cell_prob_1 = self.cell_prob_1_50d
            return self.cell_prob_1_50d
        else:
            self.cell_prob_2 = self.cell_prob_2_50d
            return self.cell_prob_2_50d

    def Prob_Reads_Given_Genotype_prob(self, gt_flag):
        if gt_flag == 0:
            return self.cell_prob_0
        elif gt_flag == 1:
            return self.cell_prob_1
        else:
            return self.cell_prob_2
