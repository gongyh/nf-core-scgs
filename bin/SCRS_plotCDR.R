#!/usr/bin/env Rscript

library(ggpubr)
library(rstatix)
library(stringr)
## draw cdr

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_cdr.R <input_cdr> <output_pdf>", call.=FALSE)
}

input_cdr <- args[1] # SCRS cdr data
output_pdf <- args[2] # store cdr ploting

pdf(output_pdf, width=8, height=5, useDingbats = FALSE)
par(ps=12)

df <- read.table(input_cdr, sep='\t', header=T)
df <- df[df$SNR>2.5,]

#df$batch <- ifelse(str_detect(df$ID_Cell, "D2O"),"batch2","batch1")

ggdotplot(df, x = "Group", y = "CDR", add = c("violin","mean_sd"), size=0.1, color="batch") +
    theme_bw() + xlab("Group") + ylab("CDR") + stat_compare_means() +
    theme(text=element_text(size=12), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

summary.stats <- df %>% select(c("Group","CDR")) %>% group_by(Group) %>% get_summary_stats(type = "common")

ggsummarytable( summary.stats, x = "Group", y = c("n", "min", "max", "mean", "median", "iqr", "sd", "se", "ci"), 
       digits=3, size=6, ggtheme = theme_bw()+
       theme(text=element_text(size=12), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
       axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) )

dev.off()

