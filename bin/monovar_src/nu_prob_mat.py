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


class Prob_matrix:

    def __init__(self, sngl_cell_ftr_list):
        self.n_cells = len(sngl_cell_ftr_list)
        self.denom_prob_matrix = np.zeros((2 * self.n_cells + 1, self.n_cells))
        self.sngl_cell_ftr_list = sngl_cell_ftr_list


    def fill_matrix(self, sngl_cell_ftr_list, original_n_cells, max_depth,
                nCr_matrix, pad):
        # Base cases (between Eq. 9 and Eq. 10)
        self.denom_prob_matrix[0, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(0, max_depth, pad))
        self.denom_prob_matrix[2, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(2, max_depth, pad))        
        self.denom_prob_matrix[1, 0] = 2 * float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(1, max_depth, pad))

        # Eq. 9
        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype(0, max_depth, pad)
            cell_j_prob_2 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype(2, max_depth, pad)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype(1, max_depth, pad)
            for l in range(0, 2 * self.n_cells + 1):
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if (l == 0):
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif (l == 1):
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    # if (sngl_cell_ftr_list[j].depth == 0):
                    # 	self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j-1]
                    # else:
                    self.denom_prob_matrix[l, j] = t1 * cell_j_prob_2 + \
                        t2 * cell_j_prob_1 + \
                        t3 * cell_j_prob_0

        for l in range(0, 2 * self.n_cells + 1):
            # Eq. 10
            self.denom_prob_matrix[l, self.n_cells - 1] = \
                self.denom_prob_matrix[l, self.n_cells - 1] \
                    / nCr_matrix[2 * self.n_cells, l]

        return self.denom_prob_matrix


    def fill_matrix_50d(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        self.denom_prob_matrix[0, 0] = float(sngl_cell_ftr_list[0] \
            .Prob_Reads_Given_Genotype_50d(0))
        self.denom_prob_matrix[2, 0] = float(sngl_cell_ftr_list[0] \
            .Prob_Reads_Given_Genotype_50d(2))
        self.denom_prob_matrix[1, 0] = 2 * float(sngl_cell_ftr_list[0] \
            .Prob_Reads_Given_Genotype_50d(1))

        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_50d(0)
            cell_j_prob_2 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_50d(2)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_50d(1)
            for l in range(0, 2 * self.n_cells + 1):
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if (l == 0):
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif (l == 1):
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    # if (sngl_cell_ftr_list[j].depth == 0):
                    # 	self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j-1]
                    # else:
                    self.denom_prob_matrix[l, j] = t1 * cell_j_prob_2 + \
                        t2 * cell_j_prob_1 + \
                        t3 * cell_j_prob_0

        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] = self.denom_prob_matrix[
                l, self.n_cells - 1] / nCr_matrix[2 * self.n_cells, l]

        return self.denom_prob_matrix


    def fill_matrix_stable(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        self.denom_prob_matrix[0, 0] = sngl_cell_ftr_list[0] \
            .Prob_Reads_Given_Genotype_prob(0)
        self.denom_prob_matrix[2, 0] = sngl_cell_ftr_list[0] \
            .Prob_Reads_Given_Genotype_prob(2)
        self.denom_prob_matrix[1, 0] = \
            2 * sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_prob(1)
        sum_l = self.denom_prob_matrix[0, 0] + self.denom_prob_matrix[2, 0] \
            + self.denom_prob_matrix[1, 0]

        self.denom_prob_matrix[0, 0] = self.denom_prob_matrix[0, 0] / sum_l
        self.denom_prob_matrix[2, 0] = self.denom_prob_matrix[2, 0] / sum_l
        self.denom_prob_matrix[1, 0] = self.denom_prob_matrix[1, 0] / sum_l

        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_prob(0)
            cell_j_prob_2 = sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_prob(2)
            cell_j_prob_1 = \
                2 * sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_prob(1)
            sum_l = 0
            for l in range(0, 2 * self.n_cells + 1):

                if (l > 2 * (j + 1)):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if (l == 0):
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif (l == 1):
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    # if (sngl_cell_ftr_list[j].depth == 0):
                    # 	self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j-1]
                    # else:
                    self.denom_prob_matrix[l, j] = t1 * cell_j_prob_2 + \
                        t2 * cell_j_prob_1 + \
                        t3 * cell_j_prob_0
                    sum_l += self.denom_prob_matrix[l, j]
            for l in range(0, 2 * (j + 1)):
                self.denom_prob_matrix[
                    l, j] = self.denom_prob_matrix[l, j] / sum_l

        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] = self.denom_prob_matrix[
                l, self.n_cells - 1] / nCr_matrix[2 * self.n_cells, l]

        return self.denom_prob_matrix
