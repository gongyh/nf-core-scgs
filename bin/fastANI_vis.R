#!/usr/bin/env Rscript

# Visualize fastANI one to one genome comparison.
# This script was modified from https://github.com/ParBLiSS/FastANI/blob/master/scripts/visualize.R

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: fastANI_vis.R <query.fasta> <subject.fasta> <visualization.visual>", call.=FALSE)
}

#Parse command line arguments
query_fasta <- args[1]
subject_fasta <- args[2]
fastANI_visual_file <- args[3]

if (!require("genoPlotR")) {
  install.packages("genoPlotR", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("genoPlotR")
}

#Read fastANI output
comparison <- try(read_comparison_from_blast(fastANI_visual_file))

#Read sequences into genoPlotR objects
Query <- try(read_dna_seg_from_file(query_fasta))
Ref <- try(read_dna_seg_from_file(subject_fasta))

plotTitle = paste(basename(query_fasta), basename(subject_fasta), sep=" v/s ")

pdf(paste0(fastANI_visual_file,".pdf"), useDingbats = FALSE, width = 10, height = 6)

plot_gene_map(dna_segs=list(Query, Ref), comparisons=list(comparison), main=plotTitle, scale=FALSE, scale_cex=1, n_scale_ticks=4)

dev.off()
