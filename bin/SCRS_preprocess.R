#!/usr/bin/env Rscript

###################################################################
### Please reformat SCRS data to txt format using labspec 6########
###################################################################

if (!require("tools")) {
   install.packages("tools", dependencies=TRUE, repos='http://cloud.r-project.org/')
   library("tools")
}

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: SCRS_preprocess.R <input_folder> <output_dir> <metadata> {<normalize_500-2000>}", call.=FALSE)
}

filepath <- file_path_as_absolute(args[1]) # directory containing SCRS txts
output <- args[2] # directory to store outputs
meta_fp <- file_path_as_absolute(args[3]) # meta TSV file, ID_Cell	Timepoint	Label	CellBg	Cell
txt_filename <- 'SCRS'

normalize_fp <- FALSE
if (length(args) >= 4) {
  normalize_fp <- as.logical(args[4])
}

dir.create(output)
output <- file_path_as_absolute(output)

# Load / install packages
if (!require("hyperSpec")) {
  install.packages("hyperSpec", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("hyperSpec")
}

if (!require("baseline")) {
  install.packages("baseline", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("baseline")
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("RColorBrewer")
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
# Read meta data          #
###########################
meta_data <- read.table(meta_fp, header=T, sep="\t", stringsAsFactors = F) # metadata table
ID_Cells <- meta_data$ID_Cell # Cells
meta_colnames <- colnames(meta_data)

###########################
# Combine txts to one txt #
###########################
setwd(filepath)
files <- list.files(pattern="*.txt") # at least 1 file
shift <- read.table(files[1], header=F, sep="\t")$V1
final_data <- c(meta_colnames, shift)
for (filename in files)
{
  ID_Cell <- sub(".txt","",filename)
  if (ID_Cell %in% ID_Cells) {
    dt <- read.table(filename, header=F, sep="\t")$V2
    # remove Cosmic Rays
    dt2 <- removeCosmic(dt)
    if (dt2$cosmic) {cat(filename,"contains cosmic ray signal.\n")}
    data <- c(unlist(meta_data[meta_data$ID_Cell==ID_Cell,]), dt2$spc)
    final_data <- cbind(final_data, data)
  }
}
setwd(output)
write.table(file=paste0(txt_filename,"_rawdata.txt"),t(final_data), sep="\t", 
            quote=F, row.names=F, col.names = F)

## transform to data frame
raw.data <- read.table(paste0(txt_filename,"_rawdata.txt"),header = T,sep="\t")

#########################
##wholespectra data######
#########################
Group <- paste(raw.data$Timepoint, sep="")
raw.data <- cbind(No_Cell=seq(1:nrow(raw.data)),Group=Group,raw.data) #Add No_Cell, Group and other meta information

### delete 1049 peak ####
#raw.data<-raw.data[,-(which(colnames(raw.data)=="X1026.09"):which(colnames(raw.data)=="X1075.76"))]

#######
##-bg##
#######
ncol_meta <- which(colnames(raw.data)=="Cell")#cols of meta data
ncol_raw.data <- ncol(raw.data)#cols of raw.data

#remove background
Cells_bgsub <- raw.data[which(raw.data$CellBg=="Cell"),]
#Cells_bgsub<-raw.data_bgsub[which(raw.data_bgsub$CellBg=="Cell"),]
#write.csv(Cells_bgsub,paste(output,"Cells_bg.csv",sep=""),quote = F,row.names = F)

#### perpare hyperSpec object ####
wavelength <- shift
data_hyperSpec <- new ("hyperSpec", data=data.frame(Cells_bgsub[,1:ncol_meta]),
                     spc = Cells_bgsub[,(ncol_meta+1):ncol_raw.data], wavelength=wavelength)

### smooth ###
data_hyperSpec <- spc.loess(data_hyperSpec,wavelength, normalize=F)

############
##baseline##
############
#data_baseline <- data_hyperSpec-spc.fit.poly.below(data_hyperSpec, data_hyperSpec, poly.order = 7)
b_als <- baseline(data_hyperSpec$spc, method='als')
data_baseline <- new ("hyperSpec", data=data.frame (Cells_bgsub[,1:ncol_meta]),
                      spc=getCorrected(b_als), wavelength=wavelength)
write.csv(data_baseline,"Cells_bg_baseline.csv",quote = F,row.names = F)

## Replace negative intensities to zero ##
data_baseline_zero <- data_baseline$spc
data_baseline_zero[data_baseline_zero<0] <- 0
data_baseline_zero_hyperSpec <- new ("hyperSpec", data=data.frame (Cells_bgsub[,1:ncol_meta]),
                                     spc = data_baseline_zero, wavelength=wavelength)
write.csv(data_baseline_zero_hyperSpec,"Cells_bg_baseline_zero.csv", quote=F, row.names=F)

#output txts
Cells_bg_baseline_zero <- "Cells_bg_baseline_zero/"
dir.create(Cells_bg_baseline_zero)
for (i in (1:nrow(data_baseline_zero_hyperSpec))){
  Cells <- data.frame(shift=shift,intensity=t(data_baseline_zero_hyperSpec[i]$spc))
  write.table(Cells,paste0(Cells_bg_baseline_zero,data_baseline_zero_hyperSpec$ID_Cell[i],".txt"),
              row.names=F,col.names=F,quote=F,sep = "\t")
}

#################
##normalization##
#################
if (normalize_fp) {
  data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowMeans (data_baseline_zero_hyperSpec[,,c(500~2000)])
} else {
  data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowMeans (data_baseline_zero_hyperSpec)
}

### remove spectrum which contains cartenoid peaks ###
keep <- c()
for (i in (1:nrow(data_baseline_zero_scale_hyperSpec))){
    peaks <- wavelength[findPeaks(data_baseline_zero_scale_hyperSpec[i]$spc)]
    if ((length(peaks[(peaks>992)&(peaks<1010)])>=1) && (length(peaks[(peaks>1145)&(peaks<1160)])>=1) &&
        (length(peaks[(peaks>1480)&(peaks<1525)])>=1)) { # cartenoid peaks
        cat(as.character(data_baseline_zero_scale_hyperSpec$ID_Cell)[i],"contains cartenoid peaks.\n")
    } else {
        keep <- c(keep, i)
    }
}
data_baseline_zero_scale_hyperSpec <- data_baseline_zero_scale_hyperSpec[keep]

#data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowSums (data_baseline_zero_hyperSpec)
write.csv(data_baseline_zero_scale_hyperSpec, "Cells_bg_baseline_zero_scale.csv",
          quote=F, row.names=F)

#output txts
Cells_bg_baseline_zero_scale <- "Cells_bg_baseline_zero_scale/"
dir.create(Cells_bg_baseline_zero_scale)
for (i in (1:nrow(data_baseline_zero_scale_hyperSpec))){
  Cells <- data.frame(shift=shift,intensity=t(data_baseline_zero_scale_hyperSpec[i]$spc))
  write.table(Cells,paste0(Cells_bg_baseline_zero_scale,data_baseline_zero_scale_hyperSpec$ID_Cell[i],".txt"),
              row.names=F,col.names=F,quote=F,sep = "\t")
}

######## calc SNR and filter by SNR ##########
SNR2 <- data_baseline_zero_scale_hyperSpec$spc
wave_nums <- shift
SNR_All <- NULL
for (i in (1:nrow(SNR2)))
{
    Baseline_start<-which.min(abs(wave_nums-1730))#1760
    Baseline_end<-which.min(abs(wave_nums-1800))#1960
    Baseline<-SNR2[i,Baseline_start:Baseline_end]
    marker<-max(SNR2[i,which.min(abs(wave_nums-3050)):which.min(abs(wave_nums-2800))]) # C-H peak
    SNR<-(marker-sum(Baseline)/length(Baseline))/sqrt(marker)
    SNR_All<-rbind(SNR_All,SNR) 
}
#str(data.frame(Cells_bgsub[,1:ncol_meta]))
#str(data.frame(SNR=SNR_All))
Meta_All <- cbind(data.frame(Cells_bgsub[keep,1:ncol_meta]), data.frame(SNR=SNR_All))

data_baseline_zero_scale_hyperSpec2 <- new ("hyperSpec", data=Meta_All, spc = SNR2, wavelength=wave_nums)

data_baseline_zero_scale_hyperSpec <- data_baseline_zero_scale_hyperSpec2[data_baseline_zero_scale_hyperSpec2$SNR>2.5]

write.csv(data_baseline_zero_scale_hyperSpec, "Cells_bg_baseline_zero_scale_hq.csv", quote=F, row.names=F)
Cells_bg_baseline_zero_scale_hq <- "Cells_bg_baseline_zero_scale_hq/"
dir.create(Cells_bg_baseline_zero_scale_hq)
for (i in (1:nrow(data_baseline_zero_scale_hyperSpec))){
   Cells <- data.frame(shift=shift,intensity=t(data_baseline_zero_scale_hyperSpec[i]$spc))
   write.table(Cells,paste0(Cells_bg_baseline_zero_scale_hq,data_baseline_zero_scale_hyperSpec$ID_Cell[i],".txt"),
               row.names=F,col.names=F,quote=F,sep = "\t")
}

#############
##   mean  ##
#############
#mean_pm_sd returns a vector with 3 values: mean - 1 sd, mean, mean + 1 sd
cluster_means_Group <- aggregate (data_baseline_zero_scale_hyperSpec, 
                                  data_baseline_zero_scale_hyperSpec$Group, mean)
write.csv(cluster_means_Group,"Cells_bg_baseline_zero_scale_M_Group.csv",quote = F,row.names = F)

#output txts
Cells_bg_baseline_zero_scale_M_Group="Cells_bg_baseline_zero_scale_M_Group/"
dir.create(Cells_bg_baseline_zero_scale_M_Group)
for (i in (1:nrow(cluster_means_Group))){
  Cells<-data.frame(shift=shift,intensity=t(cluster_means_Group[i]$spc))
  write.table(Cells,paste(Cells_bg_baseline_zero_scale_M_Group,cluster_means_Group$Group[i],"_M.txt",sep=""),row.names=F,col.names=F,quote=F,sep = "\t")
}
###################################################################

#############
## mean_sd ##
#############
#mean_pm_sd returns a vector with 3 values: mean - 1 sd, mean, mean + 1 sd
cluster_meansd_Group <- aggregate (data_baseline_zero_scale_hyperSpec, 
                                   data_baseline_zero_scale_hyperSpec$Group, mean_pm_sd)
write.csv(cluster_meansd_Group,"Cells_bg_baseline_zero_scale_Msd_Group.csv",quote = F,row.names = F)
###################################################################

#----------------------------------------------------------------------------------------------------
# output plots
#----------------------------------------------------------------------------------------------------
pdf("Meansd_spectra.pdf", width = 5,height = 5)
cluster.cols<-brewer.pal( 6, 'Dark2')
plot (cluster_meansd_Group,
      stacked = ".aggregate",
      fill = ".aggregate",
      title.args = list (xlab = "Raman shift (cm-1)",
                         col.lab = "black"),
      col = cluster.cols)
#abline(v = c(1006.48, 1159.15, 1522.86), lty=3,col="grey50")
dev.off()

pdf("Mean_spectra.pdf", width = 5,height = 5)
plot (cluster_means_Group,
      stacked = ".aggregate",
      fill = ".aggregate",
      title.args = list (xlab = "Raman shift (cm-1)",
                         col.lab = "black"),
      col = cluster.cols)
#abline(v = c(1006.48, 1159.15, 1522.86), lty=3,col="grey50")
dev.off()
#----------------------------------------------------------------------------------------------------
