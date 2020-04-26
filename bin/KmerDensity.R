#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: KmerDensity.R <working_dir> <prefix>", call.=FALSE)
}

wd <- args[1]
prefix <- args[2]

setwd(wd)

if (!require("magicaxis")) {
  install.packages("magicaxis", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("magicaxis")
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("RColorBrewer")
}

df <- data.frame()
for (i in 1:10) {
  fname <- paste0(prefix,"_cov31_p",i,".csv")
  data <- read.csv(fname, header=T)
  data$group <- paste0(i*10, "%")
  df <- rbind(df, data)
}

cols <- brewer.pal(n = 10, name = "Set3")
pNK <- c(0)
pdf(paste0(prefix,"_kmer.pdf"), useDingbats=FALSE, width=6, height=6)

opar <- par()

par(ps=10, mgp=c(1.5,1,0))

data <- df[df$group=="10%",]
maxNK <- max(df$NumKmers)
lgds <- c("10%")
pNK <- c(pNK, data[which.max(data$NumKmers[2:100]),"NumKmers"])
plot(data$Covg, data$NumKmers, main="", xlab="Kmer coverage", ylab="Number of Kmers", lwd=2,
     xlim=c(0,100), ylim=c(0, maxNK*1.05),axes=FALSE, col=cols[1],xaxs="i",yaxs="i",type="l")

magaxis(1:4)

for (i in 2:10) {
  group <- paste0(i*10, "%")
  lgds <- c(lgds, group)
  data <- df[df$group==group,]
  lines(data$Covg,data$NumKmers, col=cols[i],new=F,lty=1, lwd=2)
  pNK <- c(pNK, data[which.max(data$NumKmers[2:100]),"NumKmers"])
}

legend(x=75, y=maxNK*0.9, legend=lgds, col=cols, lty=1, lwd=2, box.col="white")

par(opar)
plot((0:10)*10, pNK, type='b', xlab="Fraction of reads (%)", ylab="# Kmers at Peak Coverage")

dev.off()

