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

Aneufinder(inputfolder=input_fn, outputfolder=output_fn, numCPU=1, binsizes=1000, pairedEndReads=TRUE,
  remove.duplicate.reads=TRUE, min.mapq=20, method='edivisive')

files <- list.files(paste0(output_fn,"/MODELS/method-edivisive/"), full.names=TRUE)
cl <- clusterByQuality(files, measures=c('spikiness','num.segments','entropy','bhattacharyya','sos'))
plot(cl$Mclust, what='classification')

selected.files <- unlist(cl$classification)
heatmapGenomewide(selected.files)

