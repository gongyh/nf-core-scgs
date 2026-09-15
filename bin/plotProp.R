#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plotProp.R <input_per-base_coverage> <prefix>", call.=FALSE)
}

input_fn <- args[1]
output_prefix <- args[2]

# Load / install packages
if (!require("magicaxis")) {
  install.packages("magicaxis", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("magicaxis")
}

data <- read.table(input_fn, header=F)$V1
data <- data[data > 1e-10]
d <- density(data)

pdf(paste0(output_prefix,"_pdrc.pdf"), useDingbats = FALSE, width = 6, height = 6)

plot(d, main=paste0("PDRC-",output_prefix), xlab="Relative Coverage", ylab="Probability Density",
    lwd=2, log="x", axes=FALSE, col="red",xaxs="i",yaxs="i")
magaxis(1:4)

dev.off()
