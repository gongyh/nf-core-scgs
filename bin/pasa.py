#!/usr/bin/env python

import logging
import traceback
import argparse
from pangraph import PanGraph

# Create an argument parser
parser = argparse.ArgumentParser(description="Run Pasa pipeline to scaffold based on pangenomes.")

# Add arguments with default values
parser.add_argument("--data_dir", required=True, help="Directory for the data")
parser.add_argument("--incomplete_sample_name", required=True, help="Name for the incomplete sample")
parser.add_argument("--assem_dir", required=True, help="Directory for the assembly")
parser.add_argument(
    "--fasta_gen", default="all", choices=["all", "partial"], help="Type of FASTA generation (default: all)"
)
parser.add_argument(
    "--PASAversion",
    default="sensitive",
    choices=["sensitive", "strict"],
    help="PASA version to use (default: sensitive)",
)
parser.add_argument("--MLR", type=int, default=0, help="PASA MLR parameter (default: 0)")
parser.add_argument("--SInfer", type=int, default=0, help="PASA SInfer parameter (default: 0)")
parser.add_argument("--min_weight_val", type=float, default=0.3, help="Minimum weight value for PASA (default: 0.3)")
parser.add_argument("--output_fasta", required=True, help="Genome sequences after scaffolding")

# Parse the arguments
args = parser.parse_args()

# Use the parsed arguments in the script
data_dir = args.data_dir
incomplete_sample_name = args.incomplete_sample_name
assem_dir = args.assem_dir
fasta_gen = args.fasta_gen
PASAversion = args.PASAversion
MLR = args.MLR
SInfer = args.SInfer
min_weight_val = args.min_weight_val
pangraph_output_sensitive = args.output_fasta

from multiprocessing import freeze_support

if __name__ == "__main__":

    freeze_support()

    pangraph = PanGraph(sample_info=None, gene_info=None, gene_position=None)

    try:
        maximum_matching = "greedy"
        pangraph.run_pangraph_pipeline(
            data_dir,
            incomplete_sample_name,
            assem_dir,
            fasta_gen,
            pangraph_output_sensitive,
            maximum_matching,
            MLR,
            SInfer,
            min_weight_val,
        )
    except Exception as e:
        logging.error(traceback.format_exc())
        # Logs the error appropriately.
