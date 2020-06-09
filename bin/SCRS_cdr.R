#!/usr/bin/env Rscript

## calc C-D ratio of SCRS (after preprocess)

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_cdr.R <input_csv> <output_dir>", call.=FALSE)
}

input_csv <- args[1] # SCRS data
output_dir <- args[2] # store cdr results

SNR1<-read.table(input_csv, header = T, sep=",")
SNR2<-SNR1[,12:length(SNR1[1,])-1]
colnames(SNR2)<-formatC(as.numeric(gsub("spc.", "", colnames(SNR2))), digits=1, format="f")
wave_nums<-as.numeric(gsub("[A-Z]", "", colnames(SNR2)))
CDR_All<-NULL
for (i in (1:nrow(SNR2)))
{
  CD_start <- which.min(abs(wave_nums-2050))
  CD_end <- which.min(abs(wave_nums-2300))
  CH_start <- which.min(abs(wave_nums-3050))
  CH_end <- which.min(abs(wave_nums-2800))
  CD <- SNR2[i,CD_start:CD_end]
  CH <- SNR2[i,CH_start:CH_end]
  CDR <- sum(CD)/(sum(CD)+sum(CH))
  CDR_All<-rbind(CDR_All, CDR)
}
CDR_All<-cbind(SNR1[,1:7], data.frame(CDR=CDR_All))
write.table(CDR_All, file=paste0(output_dir,"/","CDR.txt"), sep="\t", quote=F, row.names=F)

