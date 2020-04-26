#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plotPreSeq.R <output_of_preseq> <prefix>", call.=FALSE)
}

input_fn <- args[1]
output_prefix <- args[2]

data <- read.table(input_fn, header=T)
x <- data[,1]
y <- data[,2]
xl <- colnames(data)[1]
yl <- colnames(data)[2]

pdf(paste0(output_prefix,".pdf"), useDingbats = FALSE, width = 6, height = 6)

plot(x,y, xlab=xl, ylab=yl, lwd=1, col="red", type="b")

dev.off()

