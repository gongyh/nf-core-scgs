#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: aneuf.r <input_folder> <output_dir>", call.=FALSE)
}
input_fn <- args[1]
output_fn <- args[2]

# Load / install packages
if (!require("AneuFinder")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("AneuFinder")
  library("AneuFinder")
}

Aneufinder(inputfolder=input_fn, outputfolder=output_fn,
  numCPU=3, method='edivisive')

