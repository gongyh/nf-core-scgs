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

class Genotype_Prob_matrix:

    def __init__(self, sngl_cell_ftr_list):
        self.sngl_cell_ftr_list = sngl_cell_ftr_list
        self.n_cells = len(sngl_cell_ftr_list)
        self.denom_prob_matrix = np.zeros((2 * self.n_cells + 1, self.n_cells))


    def fill_matrix(self, nCr_matrix):
        self.denom_prob_matrix[0, 0] = float(
            self.sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_prob(0))
        self.denom_prob_matrix[2, 0] = float(
            self.sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_prob(2))
        self.denom_prob_matrix[1, 0] = \
            2 * float(self.sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_prob(1))

        for j in range(1, self.n_cells):
            cell_j_prob_0 = self.sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_prob(0)
            cell_j_prob_2 = self.sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_prob(2)
            cell_j_prob_1 = 2 * self.sngl_cell_ftr_list[j] \
                .Prob_Reads_Given_Genotype_prob(1)
            for l in range(0, 2 * self.n_cells + 1):
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if l == 0:
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif l == 1:
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]

                    self.denom_prob_matrix[l, j] = t1 * cell_j_prob_2 + \
                        t2 * cell_j_prob_1 + \
                        t3 * cell_j_prob_0

        for l in range(2 * self.n_cells + 1):
            self.denom_prob_matrix[l, -1] /= nCr_matrix[2 * self.n_cells, l]

        return None