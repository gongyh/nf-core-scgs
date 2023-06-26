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


class Single_cell_genotype_records:

    def __init__(self, current_cell_ftr_list, other_cells_ftr_list,
                curr_cell_residual_mat, n_cells, prior_variant_allele):
        self.current_cell_ftr_list = current_cell_ftr_list
        self.other_cells_ftr_list = other_cells_ftr_list
        self.residual_mat = curr_cell_residual_mat
        self.other_cell_list_len = len(other_cells_ftr_list)

        self.n_cells = n_cells
        self.prior_variant_allele = prior_variant_allele


    def find_coeff(self, n_cells, l, j, nCr_matrix):
        if j > l:
            return 0
        else:
            num = nCr_matrix[l, j] * nCr_matrix[2 * n_cells - l, 2 - j]
            denom = nCr_matrix[2 * n_cells, 2]
            return float(num) / denom


    def other_cells_genotype_prob(self, gt_flag, nCr_matrix):
        if self.other_cell_list_len == 0:
            return 1
        else:
            prob = 0.0
            for l in range(self.residual_mat.denom_prob_matrix.shape[0]):
                coeff = self.find_coeff(self.n_cells, l, gt_flag, nCr_matrix)
                prob += coeff \
                    * self.residual_mat.denom_prob_matrix[l, -1] \
                    * self.prior_variant_allele[l + gt_flag]
            return prob


    def find_genotype_prob(self, gt, nCr_matrix):
        p1 = self.current_cell_ftr_list.Prob_Reads_Given_Genotype_prob(gt)
        p2 = self.other_cells_genotype_prob(gt, nCr_matrix)
        if p1 == 0:
            p1 = 1e-322
        if p2 == 0:
            p2 = 1e-322
        return p1 * p2


    def get_genotype_prob(self, nCr_matrix):
        prob = np.zeros(3)
        for gt in (0, 1, 2):
            prob[gt] = self.find_genotype_prob(gt, nCr_matrix)
        return prob
