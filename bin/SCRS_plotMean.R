#!/usr/bin/env Rscript

if (!require("tools")) {
   install.packages("tools", dependencies=TRUE, repos='http://cloud.r-project.org/')
   library("tools")
}

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_plotMean.R <input_csv> <output_dir>", call.=FALSE)
}

input_csv <- file_path_as_absolute(args[1]) # Cells_bg_baseline_zero_scale.csv
output <- args[2] # directory to store outputs

dir.create(output)
output <- file_path_as_absolute(output)

# Load / install packages
if (!require("hyperSpec")) {
  install.packages("hyperSpec", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("hyperSpec")
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("RColorBrewer")
}

setwd(output)

### Read spectra from input csv ###
rawdata <- read.csv(input_csv, header=T)
spc <- rawdata[,9:length(rawdata[1,])-1]
colnames(spc)<-formatC(as.numeric(gsub("spc.", "", colnames(spc))), digits=1, format="f")
wave_length<-as.numeric(gsub("[A-Z]", "", colnames(spc)))
data_baseline_zero_scale_hyperSpec <- new ("hyperSpec", data=rawdata[,1:7],
                                      spc = spc, wavelength=wave_length)

######## calc SNR ##########
SNR2 <- spc
wave_nums <- wave_length
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
SNR_All<-cbind(rawdata[,1:7], data.frame(SNR=SNR_All))

data_baseline_zero_scale_hyperSpec2 <- new ("hyperSpec", data=SNR_All,
                                      spc = spc, wavelength=wave_length)

data_baseline_zero_scale_hyperSpec <- data_baseline_zero_scale_hyperSpec2[data_baseline_zero_scale_hyperSpec2$SNR>2.5]

#############
#mean_pm_sd returns a vector with 3 values: mean - 1 sd, mean, mean + 1 sd
cluster_means_Group <- aggregate (data_baseline_zero_scale_hyperSpec, 
                                  data_baseline_zero_scale_hyperSpec$Group, mean)
#write.csv(cluster_means_Group,"Cells_bg_baseline_zero_scale_M_Group.csv",quote = F,row.names = F)
###################################################################

#############
## mean_sd ##
#############
#mean_pm_sd returns a vector with 3 values: mean - 1 sd, mean, mean + 1 sd
cluster_meansd_Group <- aggregate (data_baseline_zero_scale_hyperSpec, 
                                   data_baseline_zero_scale_hyperSpec$Group, mean_pm_sd)
#write.csv(cluster_meansd_Group,"Cells_bg_baseline_zero_scale_Msd_Group.csv",quote = F,row.names = F)
###################################################################

#----------------------------------------------------------------------------------------------------
# output plots
#----------------------------------------------------------------------------------------------------
pdf("Meansd_spectra.pdf", width = 7, height = 5)
cluster.cols<-brewer.pal( 6, 'Dark2')
plot (cluster_meansd_Group,
            stacked = ".aggregate",
            fill = ".aggregate",
            title.args = list (xlab = "Raman shift (cm-1)",
                               col.lab = "black"),
            col = cluster.cols,
            axis.args=list(y=list(las=1)))
#abline(v = c(1006.48, 1159.15, 1522.86), lty=3,col="grey50")
dev.off()

pdf("Mean_spectra.pdf", width = 7, height = 5)
plot (cluster_means_Group,
            stacked = ".aggregate",
            fill = ".aggregate",
            title.args = list (xlab = "Raman shift (cm-1)",
                               col.lab = "black"),
            col = cluster.cols,
            axis.args=list(y=list(las=1)))
#abline(v = c(1006.48, 1159.15, 1522.86), lty=3,col="grey50")
dev.off()
#----------------------------------------------------------------------------------------------------

