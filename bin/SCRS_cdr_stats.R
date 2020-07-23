#!/usr/bin/env Rscript

# install necessary libraries
p <- c("optparse","tools","RColorBrewer","permute","ggplot2", "reshape2", "dplyr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

#--------------------------------------------------
getPermuteMatrix <- function(perm, N,  strata = NULL)
{
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
        perm <- how(nperm = perm) 
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm)))
            setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
        perm <- shuffleSet(N, control = perm)
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
        attr(perm, "control") <-
            structure(list(within=list(type="supplied matrix"),
                           nperm = nrow(perm)), class = "how")
    perm
}

#--------------------------------------------------
cumcv <- function(x) sapply(seq_along(x), function(k, z) sd(z[1:k])/mean(z[1:k]) , z = x)
cummedian <- function(x) sapply(seq_along(x), function(k, z) median(z[1:k]), z = x)
cummean <- function(x) cumsum(x) / seq_along(x)

#--------------------------------------------------
Stat_accum <- function (x, permutations = 100, stat="cumcv", raw = FALSE, collector = FALSE, subset, ...) {

    if (!missing(subset)) 
        x <- subset(x, subset)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    pmat <- getPermuteMatrix(permutations, n)
    permutations <- nrow(pmat)
    result <- array(dim = c(n-1, permutations))
    dimnames(result) <- list(cells = c(2:n), permutation = c(1:permutations))
    for (k in 1:permutations) {
        if(stat=="cumcv"){
        result[, k] <- cumcv(x[pmat[k, ]])[-1] 
        }else if (stat=="cummean"){
        result[, k] <- cummean(x[pmat[k, ]])[-1]
        }else if (stat=="cummax"){
        result[, k] <- cummax(x[pmat[k, ]])[-1]
        }
    }
    if (raw) 
        collector <- FALSE
    if (collector) 
        ref <- cumcv(x)[-1] 
    if (raw) {
            result <- result
    } else {
        tmp <- array(dim = c(n-1, 7 + as.numeric(collector)))
        for (i in 1:(n-1)) {
                tmp[i, 1] <- mean(result[i, 1:permutations])
                tmp[i, 2] <- median(result[i, 1:permutations])
                tmp[i, 3] <- sd(result[i, 1:permutations])
                tmp[i, 4] <- min(result[i, 1:permutations])
                tmp[i, 5] <- max(result[i, 1:permutations])
                tmp[i, 6] <- quantile(result[i, 1:permutations], 0.025)
                tmp[i, 7] <- quantile(result[i, 1:permutations], 0.975)
             if (collector) 
                tmp[i, 8] <- ref[i]
           
        }
        result <- tmp
        dimnames(result) <- list(cells = c(2:n),  
            c("mean", "median", "stdev", "min", "max", "Qnt 0.025", "Qnt 0.975", 
                if (collector) "Collector"))
    }
    #attr(result, "control") <- attr(pmat, "control")
    class(result) <- c("Stat_accum", class(result))
    result
}

#--------------------------------------------------
StatAccumCurve<-function(x, stat="cumcv", permutations = 100, Group, Group_name, outdir, width=7, height=7, scales="fixed"){
                 require("ggplot2")
                 Get_breaks<- function(Min, Max, perc){ int <- (round(Max, -(nchar(Max)-1)))*perc 
                                                       breaks <- c(Min, seq(from=int, to=Max, by=int))
                                                       return(breaks)
                       }
                 x <- as.matrix(x)
                 if (missing(Group)) {                       
                       ra <- Stat_accum(x, stat=stat, permutations=permutations)
                       ra_df <- data.frame(cells = 2:(nrow(ra)+1), ra) 
                       p <- ggplot(ra_df, aes(x = cells, y = mean)) + geom_point(size=0.1) + geom_line() +
                            xlab("Number of cells pooled") + ylab( stat ) +
                            #geom_ribbon(data = ra_df, aes(ymin=Qnt.0.025,ymax=Qnt.0.975), alpha=0.2)+
                            geom_ribbon(data = ra_df, aes(ymin=mean-stdev,ymax=mean+stdev), alpha=0.2)+
                            geom_vline(aes(xintercept=3), colour="#990000", linetype="dashed")+
                            geom_vline(aes(xintercept=20), colour="#990000", linetype="longdash")+
                            geom_vline(aes(xintercept=60), colour="#990000", linetype="solid")+
                            theme_bw()
                       ggsave(filename=paste(outdir, stat, "_accum.pdf", sep=""), plot=p, width=4, height=4)
                       sink(paste(outdir, stat, "_accum_df.txt", sep="")); write.table(ra_df, quote=FALSE, sep="\t", row.names=FALSE); sink()
                 }else{
                       x_list<-split(x, Group)
                       ra_list<-lapply(names(x_list), function(x){ ra<-Stat_accum(x_list[[x]], stat=stat, permutations=permutations); ra_df <- data.frame(group=rep(x, nrow(ra)), cells = 2:(nrow(ra)+1), ra)})
                       ra_df<-do.call(rbind, ra_list)
                       if(length(Group_name)==1 & length(Group)==nrow(x)){
                         colnames(ra_df)[1]<-Group_name
                         sink(paste(outdir, stat, "_accum_df_by_",Group_name,".txt", sep="")); write.table(ra_df, quote=FALSE, sep="\t", row.names=FALSE); sink()
                         }else{
                         groups<-data.frame(do.call(rbind, strsplit(as.character(ra_df$group), "\\."))); names(groups)<-Group_name
                         ra_df<-data.frame(groups, ra_df)
                         sink(paste(outdir, stat, "_accum_df_by_",Group_names,".txt", sep="")); write.table(ra_df, quote=FALSE, sep="\t", row.names=FALSE); sink()
                         }
                       
                       breaks<-Get_breaks(min(ra_df$cells), max(ra_df$cells), 1/6)
                       p <- ggplot(ra_df, aes(x = cells, y = mean)) + geom_point(size=0.2) + geom_line() +
                            xlab("Number of cells pooled") + ylab(stat) +
                            #geom_ribbon(data = ra_df, aes(ymin=Qnt.0.025,ymax=Qnt.0.975), alpha=0.2)+
                            geom_ribbon(data = ra_df, aes(ymin=mean-stdev,ymax=mean+stdev), alpha=0.2)+
                            geom_vline(aes(xintercept=3), colour="#990000", linetype="dashed")+
                            geom_vline(aes(xintercept=20), colour="#990000", linetype="longdash")+
                            geom_vline(aes(xintercept=60), colour="#990000", linetype="solid")+
                            theme_bw() + scale_x_continuous(breaks=breaks)
                       if(length(Group_name)==1 & length(Group)==nrow(x)){
                            p <- p + facet_wrap(~ get(Group_name), scales=scales)
                            ggsave(filename=paste(outdir, stat, "_accum_facets_by_",Group_name,".pdf", sep=""), plot=p, limitsize=FALSE, width=width, height=height)
                       }else{
                            p <- p + facet_grid( get(Group_name[1]) ~ get(Group_name[2]), scales=scales)
                            Group_names<-paste(Group_name, collapse="__")
                            ggsave(filename=paste(outdir, stat, "_accum_facets_by_",Group_names,".pdf", sep=""), plot=p, limitsize=FALSE, width=width, height=height)
                       }
                       
                 }
                 invisible(p)
}
#--------------------------------------------------
# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_cdr_stats.R <input_cdr_txt> <output_dir>", call.=FALSE)
}

input <- file_path_as_absolute(args[1]) # SCRS CDR results, CDR.txt
outpath <- args[2] # store rarefaction results

setwd("./")

width=7
height=7
scales="free"
options(warn=-1)

dir.create(outpath)
outpath <- file_path_as_absolute(outpath)

mat <- read.table(input, header = T,  sep="\t")
mat$CD_ratio <- mat$CDR
mat$Time <- mat$Timepoint
rownames(mat) <- mat$No_Cell

######negtive value to zero for cummean and cummax#########
mat1<-mat
mat1[which(mat1$CD_ratio<0),"CD_ratio"]<-0
#mat1[which(mat1$Protein<0),"Protein"]<-0
#mat1[which(mat1$TAG<0),"TAG"]<-0
#####figure parameters setting##############
#####stat="cummean"
stat="cummean"
for (pheno in (c("CD_ratio"))){
  x<-mat1[, pheno]; names(x)<-rownames(mat1)
  StatAccumCurve(x, stat=stat, permutations = 1000, Group=mat[, "Time"], Group_name="Time", 
       outdir=paste0(outpath, "/", pheno,"_"), width=width, height=height, scales=scales)
  StatAccumCurve(x, stat=stat, permutations = 1000, outdir=paste0(outpath, "/",pheno,"_allcells_"), 
       width=width, height=height, scales=scales)
}
####stat="cummax"
stat="cummax"
for (pheno in (c("CD_ratio"))){
  x<-mat1[, pheno]; names(x)<-rownames(mat1)
  StatAccumCurve(x, stat=stat, permutations = 1000, Group=mat[, "Time"], Group_name="Time", 
       outdir=paste0(outpath, "/", pheno,"_"), width=width, height=height, scales=scales)
  StatAccumCurve(x, stat=stat, permutations = 1000, outdir=paste0(outpath, "/", pheno,"_allcells_"), 
       width=width, height=height, scales=scales)
}
######negtive value to 1e-10 for cumcv#########
mat1<-mat
mat1[which(mat1$CD_ratio<0),"CD_ratio"]<-1e-10
#mat1[which(mat1$Protein<0),"Protein"]<-1e-10
#mat1[which(mat1$TAG<0),"TAG"]<-1e-10

stat="cumcv"
scales="free"
scales="fixed"
for (pheno in (c("CD_ratio"))){
  x<-mat1[, pheno]; names(x)<-rownames(mat1)
  StatAccumCurve(x, stat=stat, permutations = 1000, Group=mat[, "Time"], Group_name="Time", 
       outdir=paste0(outpath, "/", pheno,"_"), width=width, height=height, scales=scales)
  StatAccumCurve(x, stat=stat, permutations = 1000, outdir=paste0(outpath, "/", pheno,"_allcells_"), 
       width=width, height=height, scales=scales)
}
###############################################

