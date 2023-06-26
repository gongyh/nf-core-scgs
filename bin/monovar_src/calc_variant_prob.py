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
from nu_prob_mat import Prob_matrix


class Calc_Var_Prob():
    """docstring for  Calc_Var_Prob"""

    def __init__(self, read_supported_cell_list):
        self.read_supported_cell_list = read_supported_cell_list
        self.matrix = Prob_matrix(read_supported_cell_list)


    def calc_zero_var_prob(self, n_cells, max_depth, nCr_matrix, pad,
                prior_variant_number):
        self.matrix.denom_prob_matrix = self.matrix.fill_matrix(
            self.read_supported_cell_list, n_cells, max_depth, nCr_matrix, pad)

        probs = self.matrix.denom_prob_matrix[:, -1] * prior_variant_number
        # Numerator Eq. 6: P(D|l=0) * p(l=0) 
        numerator = probs[0]
        # Denominator Eq. 6: \sum_{l=0}^{2m} P(D|l) * p(l) 
        denominator = np.sum(probs)

        if denominator == 0.0:
            self.matrix.denom_prob_matrix = self.matrix.fill_matrix_50d(
                self.read_supported_cell_list, n_cells, nCr_matrix)
            probs = self.matrix.denom_prob_matrix[:, -1] * prior_variant_number
            # Numerator Eq. 6: P(D|l=0) * p(l=0) 
            numerator = probs[0]
            # Denominator Eq. 6: \sum_{l=0}^{2m} P(D|l) * p(l) 
            denominator = np.sum(probs)

            if denominator == 0:
                self.matrix.denom_prob_matrix = self.matrix.fill_matrix_stable(
                    self.read_supported_cell_list, n_cells, nCr_matrix)
                probs = self.matrix.denom_prob_matrix[:, -1] \
                    * prior_variant_number
                # Numerator Eq. 6: P(D|l=0) * p(l=0) 
                numerator = probs[0]
                # Denominator Eq. 6: \sum_{l=0}^{2m} P(D|l) * p(l) 
                denominator = np.sum(probs)

                if denominator == 0:
                    denominator = 1e-300

        zero_variant_prob = numerator / denominator
        return zero_variant_prob, denominator
