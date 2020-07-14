#!/usr/bin/env Rscript

## calc SNR of SCRS (after preprocess)

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_snr.R <input_csv> <output_dir>", call.=FALSE)
}

input_csv <- args[1] # SCRS data
output_dir <- args[2] # store snr results

dir.create(output_dir)

SNR1<-read.table(input_csv, header = T, sep=",")
SNR2<-SNR1[,12:length(SNR1[1,])-1]
colnames(SNR2)<-formatC(as.numeric(gsub("spc.", "", colnames(SNR2))), digits=1, format="f")
wave_nums<-as.numeric(gsub("[A-Z]", "", colnames(SNR2)))
SNR_All<-NULL
for (i in (1:nrow(SNR2)))
{
  Baseline_start<-which.min(abs(wave_nums-1730))#1760
  Baseline_end<-which.min(abs(wave_nums-1800))#1960
  Baseline<-SNR2[i,Baseline_start:Baseline_end]
  marker<-max(SNR2[i,which.min(abs(wave_nums-3050)):which.min(abs(wave_nums-2800))]) # C-H peak
  SNR<-(marker-sum(Baseline)/length(Baseline))/sqrt(marker)
  SNR_All<-rbind(SNR_All,SNR) 
}
SNR_All<-cbind(SNR1[,1:7], data.frame(SNR=SNR_All))
write.table(SNR_All, file=paste0(output_dir,"/","SNR.txt"), sep="\t", quote=F, row.names=F)

