#!/usr/bin/env Rscript

library(ggpubr)
library(rstatix)

## calc C-D ratio of SCRS (after preprocess)

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_cdr.R <input_csv> <output_dir>", call.=FALSE)
}

input_csv <- args[1] # SCRS data
output_dir <- args[2] # store cdr results

dir.create(output_dir)

SNR1<-read.table(input_csv, header = T, sep=",")
SNR2<-SNR1[,12:length(SNR1[1,])-1]
colnames(SNR2)<-formatC(as.numeric(gsub("spc.", "", colnames(SNR2))), digits=1, format="f")
wave_nums<-as.numeric(gsub("[A-Z]", "", colnames(SNR2)))
SNR_All<-NULL
CDR_All<-NULL
for (i in (1:nrow(SNR2)))
{
  Baseline_start<-which.min(abs(wave_nums-1730))#1760
  Baseline_end<-which.min(abs(wave_nums-1800))#1960
  Baseline<-SNR2[i,Baseline_start:Baseline_end]
  marker<-max(SNR2[i,which.min(abs(wave_nums-3050)):which.min(abs(wave_nums-2800))]) # C-H peak
  SNR<-(marker-sum(Baseline)/length(Baseline))/sqrt(marker)
  SNR_All<-rbind(SNR_All,SNR)

  CD_start <- which.min(abs(wave_nums-2050))
  CD_end <- which.min(abs(wave_nums-2300))
  CH_start <- which.min(abs(wave_nums-3050))
  CH_end <- which.min(abs(wave_nums-2800))
  CD <- SNR2[i,CD_start:CD_end]
  CH <- SNR2[i,CH_start:CH_end]
  CDR <- sum(CD)/(sum(CD)+sum(CH))
  CDR_All<-rbind(CDR_All, CDR)
}
CDR_All<-cbind(SNR1[,1:7], data.frame(SNR=SNR_All, CDR=CDR_All))
write.table(CDR_All, file=paste0(output_dir,"/","SNR_CDR.txt"), sep="\t", quote=F, row.names=F)

pdf(paste0(output_dir,"/","CDR.pdf"), width=7, height=5, useDingbats = FALSE)
par(ps=12)
df <- read.table(paste0(output_dir,"/","SNR_CDR.txt"), sep='\t', header=T)
df <- df[df$SNR>2,]

ggdotplot(df, x = "Group", y = "CDR", add = c("violin","mean_sd"), size=0.3) +
    theme_bw() + xlab("Group") + ylab("CDR") + stat_compare_means() +
    theme(text=element_text(size=12), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

summary.stats <- df %>% select(c("Group","CDR")) %>% group_by(Group) %>% get_summary_stats(type = "common")

ggsummarytable( summary.stats, x = "Group", y = c("n", "min", "max", "mean", "median", "iqr", "sd", "se", "ci"), 
       digits=3, size=6, ggtheme = theme_bw()+
       theme(text=element_text(size=12), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
       axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) )

dev.off()

