#!/usr/bin/env Rscript

###################################################################
### Please reformat SCRS data to txt format using labspec 6########
###################################################################

if (!suppressWarnings(suppressMessages(require("tools",quietly=T)))) {
   install.packages("tools", dependencies=TRUE, repos='http://cloud.r-project.org/')
   suppressWarnings(suppressMessages(library("tools")))
}

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: SCRS_preprocess.R <input_txt> <output_dir>", call.=FALSE)
}

filepath <- file_path_as_absolute(args[1]) # SCRS txt
output <- args[2] # directory to store outputs
txt_filename <- 'SCRS'

dir.create(output, showWarnings=F)
output <- file_path_as_absolute(output)

# Load / install packages
if (!suppressWarnings(suppressMessages(require("hyperSpec",quietly=T)))) {
  install.packages("hyperSpec", dependencies=TRUE, repos='http://cloud.r-project.org/')
  suppressWarnings(library("hyperSpec"))
}

if (!suppressWarnings(suppressMessages(require("baseline",quietly=T)))) {
  install.packages("baseline", dependencies=TRUE, repos='http://cloud.r-project.org/')
  suppressWarnings(library("baseline"))
}

if (!suppressWarnings(suppressMessages(require("RColorBrewer",quietly=T)))) {
  install.packages("RColorBrewer", dependencies=TRUE, repos='http://cloud.r-project.org/')
  suppressWarnings(library("RColorBrewer"))
}

sys.script <- function(){
  cmdargs <- commandArgs(trailingOnly = FALSE)
  fl <- grep("--file=", cmdargs)
  if (length(fl) > 0) {
    # Rscript
    return(normalizePath(gsub("--file=", "", cmdargs[fl])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

source(paste(dirname(sys.script()),'SCRS_helpers.R',sep="/"))

###########################
# Read SCRS txt #
###########################
df <- read.table(filepath, header=F, sep="\t")
shift <- df$V1
dt <- df$V2
# remove Cosmic Rays
dt2 <- removeCosmic(dt)
if (dt2$cosmic) {
    cat("contains cosmic ray signal.\n")
} else {
    cat("Cosmic ray signal not detected.\n")
}


setwd(output)

#### perpare hyperSpec object ####
wavelength<-shift
data_hyperSpec<-new ("hyperSpec", spc = dt2$spc, wavelength=wavelength)

### smooth ###
data_hyperSpec<-spc.loess(data_hyperSpec, wavelength, normalize=F)

############
##baseline##
############
#data_baseline <- data_hyperSpec-spc.fit.poly.below(data_hyperSpec, data_hyperSpec, poly.order = 7)
b_als <- baseline(data_hyperSpec$spc, method='als')
data_baseline <- new ("hyperSpec", spc=getCorrected(b_als), wavelength=wavelength)

## Replace negative intensities to zero ##
data_baseline_zero<-data_baseline$spc
data_baseline_zero[data_baseline_zero<0] <- 0
data_baseline_zero_hyperSpec<-new ("hyperSpec", spc = data_baseline_zero, wavelength=wavelength)

#output txts
Cells_bg_baseline_zero <- "Cells_bg_baseline_zero"
Cells<-data.frame(shift=shift,intensity=t(data_baseline_zero_hyperSpec[1]$spc))
write.table(Cells, paste0(Cells_bg_baseline_zero,".txt"), row.names=F, col.names=F, quote=F, sep = "\t")

#################
##normalization##
#################
data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowMeans (data_baseline_zero_hyperSpec)

### remove spectrum which contains cartenoid peaks ###
peaks <- wavelength[findPeaks(data_baseline_zero_scale_hyperSpec[1]$spc)]
if ((length(peaks[(peaks>992)&(peaks<1010)])>=1) && (length(peaks[(peaks>1145)&(peaks<1160)])>=1) &&
    (length(peaks[(peaks>1480)&(peaks<1525)])>=1)) { # cartenoid peaks
    cat("Contain cartenoid peaks.\n")
} else {
    cat("Cartenoid peaks not detected.\n")
}

#output txt
Cells_bg_baseline_zero_scale <- "Cells_bg_baseline_zero_scale"
Cells<-data.frame(shift=shift,intensity=t(data_baseline_zero_scale_hyperSpec[1]$spc))
write.table(Cells,paste0(Cells_bg_baseline_zero_scale,".txt"), row.names=F,col.names=F,quote=F,sep = "\t")

#### SNR and CDR ####
wave_nums <- wavelength
SNR2 <- data_baseline_zero_scale_hyperSpec[1]$spc

Baseline_start<-which.min(abs(wave_nums-1730))#1760
Baseline_end<-which.min(abs(wave_nums-1800))#1960
Baseline<-SNR2[1,Baseline_start:Baseline_end]
marker<-max(SNR2[1,which.min(abs(wave_nums-3050)):which.min(abs(wave_nums-2800))]) # C-H peak
SNR<-(marker-sum(Baseline)/length(Baseline))/sqrt(marker)
cat("SNR = ", as.character(SNR) ," \n")

CD_start <- which.min(abs(wave_nums-2050))
CD_end <- which.min(abs(wave_nums-2300))
CH_start <- which.min(abs(wave_nums-3050))
CH_end <- which.min(abs(wave_nums-2800))
CD <- SNR2[1,CD_start:CD_end]
CH <- SNR2[1,CH_start:CH_end]
CDR <- sum(CD)/(sum(CD)+sum(CH))
cat("CDR = ", as.character(CDR) ," \n")


