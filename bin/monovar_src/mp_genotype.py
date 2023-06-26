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

from genotype_prob_mat import Genotype_Prob_matrix
from nu_genotype_single_cell import Single_cell_genotype_records
import numpy as np
np.seterr(divide='ignore')

POOR_QUAL = 1e-100
PHRED_POOR_QUAL = -10 * np.log10(POOR_QUAL)
MAX_GQ = 99
genotype_dict = {0: '0/0', 1: '0/1', 2: '1/1'}


def get_info_string(read_supported_cell_list, n_cells, nCr_matrix,
            prior_variant_number, denominator, cell_count):

    current_cell_ftr_info = read_supported_cell_list[cell_count]
    cp_read_supported_cell_list = [j for i, j \
        in enumerate(read_supported_cell_list) if i != cell_count]

    curr_cell_residual_mat = Genotype_Prob_matrix(cp_read_supported_cell_list)
    if len(cp_read_supported_cell_list) > 0:
        curr_cell_residual_mat.fill_matrix(nCr_matrix)
  
    cell_genotype_obj = Single_cell_genotype_records(
        current_cell_ftr_info, cp_read_supported_cell_list, 
        curr_cell_residual_mat, n_cells, prior_variant_number
    )
    p_list = cell_genotype_obj.get_genotype_prob(nCr_matrix) / denominator
    max_p_g = p_list.max()
    # Determining GT
    if max_p_g == 0:
        p_list[0] = current_cell_ftr_info.cell_prob_0
        p_list[1] = current_cell_ftr_info.cell_prob_1
        p_list[2] = current_cell_ftr_info.cell_prob_2
        max_p_g = p_list.max()
    
    if max_p_g == 0:
        if current_cell_ftr_info.Alt_freq < 0.1:
            g_ind = 0
        elif current_cell_ftr_info.Alt_freq > 0.8:
            g_ind = 2
        else:
            g_ind = 1
        norm_p_list = np.full(3, POOR_QUAL)
        norm_p_list[g_ind] = 1
    else:
        g_ind = p_list.argmax()
        norm_p_list = p_list / p_list.sum()

    # Determining PL
    PL = np.round(-10 * np.log10(norm_p_list))
    PL = np.where(np.isinf(PL), PHRED_POOR_QUAL, PL)
    # Determining GQ
    GQ = min(MAX_GQ, np.round(-10 * np.log10(1 - norm_p_list[g_ind])))

    final_genotype = genotype_dict[g_ind]
    info = '{gt}:{ad_r},{ad_a}:{dp}:{gq:.0f}:{pl0:.0f},{pl1:.0f},{pl2:.0f}' \
        .format(gt=final_genotype, ad_r=current_cell_ftr_info.refDepth,
            ad_a=current_cell_ftr_info.alt_count,
            dp=current_cell_ftr_info.depth, gq=GQ, pl0=PL[0], pl1=PL[1],
            pl2=PL[2])
    return info, g_ind
