#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 2020

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
from scipy.stats import beta as scibeta


class doubletFinder:
    def __init__(
        self,
        df_total,
        df_alt,
        cells,
        delta,
        beta,
        mu_hom,
        alpha_fp,
        alpha_fn,
        missing=False,
        verbose=True,
        binom=False,
        precision=None,
    ):
        self.df_total = df_total
        self.df_alt = df_alt
        self.cells = cells
        self.delta = delta
        self.beta = beta
        self.mu_hom = mu_hom
        self.alpha_fp = alpha_fp
        self.alpha_fn = alpha_fn
        self.precision = precision
        self.missing = missing

        if binom:
            self.prv_y = self.prv_y_b
        else:
            self.prv_y = self.prv_y_bb

        self.muts = list(self.df_total.columns)
        npos = len(list(self.df_total.columns))
        print(f"number of position is {npos}")

        if (self.df_total.values - self.df_alt.values).min() < 0:
            raise Exception("total reads must be greater than or equal to alternate reads!")

        self.Sigma = ((0, 1), (0, 1 / 2, 1))
        self.Theta = ((0, 1), (0, 1 / 2, 1))

        self.px_z = {x: 0 for z in [0, 1] for x in itertools.product(self.muts, self.Sigma[z], [z])}

        if self.alpha_fn == None or self.alpha_fp == None or self.precision == None:
            # get vaf numpy array
            vaf_values = (self.df_alt / self.df_total).values
            vaf_values = vaf_values[~np.isnan(vaf_values)]

            # filter vafs less than 0.15 and get alpha_fp
            if self.alpha_fp == None or self.precision == None:
                mean1, prec1 = getBetaMOM(vaf_values[vaf_values <= 0.15])
                self.alpha_fp = mean1

            if self.alpha_fn == None or self.precision == None:
                mean2, prec2 = getBetaMOM(vaf_values[vaf_values >= 0.85])
                self.alpha_fn = 1 - mean2

            if self.precision == None:
                mean3, prec3 = getBetaMOM(vaf_values[(vaf_values > 0.15) & (vaf_values < 0.85)])
                self.precision = np.median([prec1, prec2, prec3])

        print(f"alpha_fn = {self.alpha_fn}")
        print(f"alpha_fp = {self.alpha_fp}")
        print(f"precision = {self.precision}")

        if self.mu_hom == None:
            estimate = True
        else:
            estimate = False

        for m in self.muts:
            if estimate:
                reads = self.df_alt[m]
                total = self.df_total[m]
                vaf = pd.DataFrame(reads / total)

                loh_cells = vaf.loc[vaf[m] >= 0.85]
                wt_cells = vaf.loc[vaf[m] <= 0.15]
                all_cells = loh_cells.shape[0] + wt_cells.shape[0]
                est_wt_rate = wt_cells.shape[0] / all_cells
                est_loh_rate = loh_cells.shape[0] / all_cells

                self.px_z[m, 0, 0] = est_wt_rate
                self.px_z[m, 1, 0] = est_loh_rate

            else:
                self.px_z[m, 0, 0] = 1 - self.mu_hom
                self.px_z[m, 1, 0] = self.mu_hom

            # norm_const = {}
            # for m in self.muts:
            #    for a, b in itertools.product(self.Sigma[0], repeat=2):
            #        c = (a + b)/2
            #        if c in self.Sigma[1]:
            #            self.px_z[m, c, 1] += self.px_z[m, a, 0] * \
            #                self.px_z[m, b, 0]
            #    norm_const[m] = sum([self.px_z[m, a, 0] * self.px_z[m, b, 0]
            #                       for a, b in itertools.product(self.Sigma[0], repeat=2)])

            # for m in self.muts:
            #     for c in self.Sigma[1]:
            #         self.px_z[m, c, 1] /= norm_const[m]

            self.px_z[m, 0, 1] = est_wt_rate * 0.25
            self.px_z[m, 1, 1] = est_loh_rate * 0.25
            self.px_z[m, 1 / 2, 1] = 0.75

        self.py_xz = {x: 0 for z in [0, 1] for x in itertools.product(self.Theta[z], self.Sigma[z], [z])}

        self.py_xz[0, 0, 0] = 1 - self.beta
        self.py_xz[1, 1, 0] = 1 - self.beta
        self.py_xz[0, 0, 1] = 1 - self.beta**2
        self.py_xz[0, 1 / 2, 1] = self.beta * (1 - self.beta)
        self.py_xz[1 / 2, 1 / 2, 1] = (1 - self.beta) ** 2
        self.py_xz[1, 1 / 2, 1] = self.beta * (1 - self.beta)
        self.py_xz[1, 1, 1] = 1 - self.beta**2

        self.doublet_result = None

        if self.delta == 0:
            self.threshold = sys.float_info.max
        elif self.delta == 1:
            self.threshold = -sys.float_info.max
        else:
            self.threshold = math.log((1 - self.delta) / self.delta)

    def solve(self):
        self.doublet_result = {}
        self.logprobs = {}
        for cell in self.cells:
            self.logprobs[cell, 0] = self.prv_z(cell, 0)
            self.logprobs[cell, 1] = self.prv_z(cell, 1)

            if self.logprobs[cell, 1] - self.logprobs[cell, 0] > self.threshold:
                self.doublet_result[cell] = "doublet"
            else:
                self.doublet_result[cell] = "singlet"

    def prv_z(self, cell, z):
        log_prob_sum = 0
        for mut in self.muts:
            v = self.df_alt.loc[cell, mut]
            r = self.df_total.loc[cell, mut] - v
            if r + v == 0:
                if z == 0 and self.beta > 0:
                    log_prob_sum += 1 * math.log(self.beta)
                if z == 1 and self.beta > 0:
                    log_prob_sum += 2 * math.log(self.beta)
                continue
            prob_sum = 0
            for x in self.Sigma[z]:
                prob_sum_x = 0
                for y in self.Theta[z]:
                    prob_sum_x += self.prv_y(r, v, y) * self.py_xz[y, x, z]
                prob_sum += self.px_z[mut, x, z] * prob_sum_x
            log_prob_sum += math.log(prob_sum)
        return log_prob_sum

    def prv_y_b(self, r, v, y):
        yprime = self.alpha_fp + (1 - self.alpha_fp - self.alpha_fn) * y
        return nCr(r + v, v) * (yprime**v) * ((1 - yprime) ** r)

    def prv_y_bb(self, r, v, y):
        yprime = self.alpha_fp + (1 - self.alpha_fp - self.alpha_fn) * y
        if y == 0:
            yprime = 0.001
        if y == 1:
            yprime = 0.999
        # print("y = " + str(y))
        alpha = self.precision * yprime
        # print("alpha = " + str(alpha))
        beta = self.precision - alpha
        # print("beta = " + str(beta))
        n = r + v
        vaf = v / n
        if vaf == 0:
            vaf = 0.01
        if vaf == 1:
            vaf = 0.99
        # print(f"n:{n} r:{r} v:{v} p:{y} alpha:{alpha} beta:{beta}")
        prob = scibeta.pdf(vaf, alpha, beta)
        return prob

    def writeSolution(self, outputFile):
        with open(outputFile, "w") as output:
            output.write("cell_id\tprob_z0\tprob_z1\tprediction\n")
            for cell in self.cells:
                output.write(
                    f"{cell}\t{self.logprobs[cell, 0]}\t{self.logprobs[cell,1]}\t{self.doublet_result[cell]}\n"
                )

    def likelihood(self):
        likelihood = len(self.cells) * math.log(1 - self.delta)
        for cell in self.cells:
            if self.doublet_result[cell] == "doublet":
                likelihood += self.logprobs[cell, 1] - self.logprobs[cell, 0] - self.threshold
        return likelihood


def getBetaMOM(x):
    m_x = np.mean(x)
    s_x = np.std(x)
    x_alpha = m_x * ((m_x * (1 - m_x) / s_x**2) - 1)
    x_beta = (1 - m_x) * ((m_x * (1 - m_x) / s_x**2) - 1)

    return x_alpha / (x_alpha + x_beta), x_alpha + x_beta


def nCr(n, r):
    f = math.factorial
    # print(n, r)
    return f(n) // f(r) // f(n - r)


def main(args):
    df_total = pd.read_csv(args.inputTotal)
    df_alt = pd.read_csv(args.inputAlternate)
    cells = list(df_total["cell_id"].values)
    df_total = df_total.set_index(["cell_id"])
    df_alt = df_alt.set_index(["cell_id"])
    vaf_values = (df_alt / df_total).values
    valid_cols = np.logical_and(np.any(vaf_values > 0.85, axis=0), np.any(vaf_values < 0.15, axis=0))
    valid_list = valid_cols.tolist()
    df_total = df_total.loc[:, valid_list]
    df_alt = df_alt.loc[:, valid_list]

    if len(df_total) != len(df_alt):
        raise Exception("number of cells in the two input files do not match!")

    if len(df_total.columns) != len(df_alt.columns):
        raise Exception("number of cells in the two input files do not match!")

    ncells = len(df_total)

    if args.verbose:
        print(f"number of cells is {ncells}")

    solver = doubletFinder(
        df_total,
        df_alt,
        cells,
        delta=args.delta,
        beta=args.beta,
        missing=args.missing,
        mu_hom=args.mu_hom,
        alpha_fp=args.alpha_fp,
        alpha_fn=args.alpha_fn,
        verbose=args.verbose,
        binom=args.binom,
        precision=args.prec,
    )

    solver.solve()

    solver.writeSolution(args.outputfile)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inputTotal", type=str, help="csv file with a table of total read counts for each position in each cell"
    )
    parser.add_argument(
        "--inputAlternate",
        type=str,
        help="csv file with a table of alternate read counts for each position in each cell",
    )
    parser.add_argument("--delta", type=float, default=0.1, help="expected doublet rate [0.1]")
    parser.add_argument("--beta", type=float, default=0.05, help="allelic dropout (ADO) rate [0.05]")
    parser.add_argument("--mu_hom", type=float, help="homozygous mutation rate [None]")
    parser.add_argument("--alpha_fp", type=float, help="copy false positive error rate [None]")
    parser.add_argument("--alpha_fn", type=float, help="copy false negative error rate [None]")
    parser.add_argument("-o", "--outputfile", type=str, help="output file name")
    parser.add_argument(
        "--noverbose",
        dest="verbose",
        help="do not output statements from internal solvers [default is false]",
        action="store_false",
    )
    parser.add_argument(
        "--binomial",
        dest="binom",
        help="use binomial distribution for read count model [default is false]",
        action="store_true",
    )
    parser.add_argument("--prec", type=float, help="precision for beta-binomial distribution [None]")
    parser.add_argument("--missing", dest="missing", help="use missing data in the model? [No]", action="store_true")
    parser.set_defaults(missing=False)
    parser.set_defaults(binom=False)
    parser.set_defaults(verbose=True)

    args = parser.parse_args(None if sys.argv[1:] else ["-h"])

    return args


def main_cli():
    args = get_options()
    main(args)


if __name__ == "__main__":
    args = get_options()
    main(args)
