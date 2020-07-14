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
if (length(args) < 4) {
  stop("Usage: SCRS_preprocess.R <input_folder> <output_dir> <metadata> <prefix>", call.=FALSE)
}

filepath <- file_path_as_absolute(args[1]) # directory containing SCRS txts
output <- args[2] # directory to store outputs
meta_fp <- file_path_as_absolute(args[3]) # meta TSV file, ID_Cell	Timepoint	Label	CellBg	Cell
txt_filename <- args[4]

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
    data <- c(unlist(meta_data[meta_data$ID_Cell==ID_Cell,]), dt)
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
Group<-paste(raw.data$Timepoint, sep="")
raw.data<-cbind(No_Cell=seq(1:nrow(raw.data)),Group=Group,raw.data) #Add No_Cell, Group and other meta information

### delete 1049 peak ####
#raw.data<-raw.data[,-(which(colnames(raw.data)=="X1026.09"):which(colnames(raw.data)=="X1075.76"))]

#######
##-bg##
#######
ncol_meta<-which(colnames(raw.data)=="Cell")#cols of meta data
ncol_raw.data<-ncol(raw.data)#cols of raw.data

#remove background
Cells_bgsub<-raw.data[which(raw.data$CellBg=="Cell"),]
#Cells_bgsub<-raw.data_bgsub[which(raw.data_bgsub$CellBg=="Cell"),]
#write.csv(Cells_bgsub,paste(output,"Cells_bg.csv",sep=""),quote = F,row.names = F)

############
##baseline##
############
wavelength<-shift
data_hyperSpec<-new ("hyperSpec", data=data.frame (Cells_bgsub[,1:ncol_meta]),
                     spc = Cells_bgsub[,(ncol_meta+1):ncol_raw.data], wavelength=wavelength)
data_baseline <- data_hyperSpec-spc.fit.poly.below(data_hyperSpec, data_hyperSpec, poly.order = 7)
#data_baseline <- data_hyperSpec - spc.rubberband(data_hyperSpec, noise=300, df=20)
write.csv(data_baseline,"Cells_bg_baseline.csv",quote = F,row.names = F)

## Replace negative intensities to zero ##
data_baseline_zero<-data_baseline$spc
data_baseline_zero[data_baseline_zero<0] <- 0
data_baseline_zero_hyperSpec<-new ("hyperSpec", data=data.frame (Cells_bgsub[,1:ncol_meta]),
                                   spc = data_baseline_zero, wavelength=wavelength)
write.csv(data_baseline_zero_hyperSpec,"Cells_bg_baseline_zero.csv", quote=F, row.names=F)

#output txts
Cells_bg_baseline_zero <- "Cells_bg_baseline_zero/"
dir.create(Cells_bg_baseline_zero)
for (i in (1:nrow(data_baseline_zero_hyperSpec))){
  Cells<-data.frame(shift=shift,intensity=t(data_baseline_zero_hyperSpec[i]$spc))
  write.table(Cells,paste0(Cells_bg_baseline_zero,data_baseline_zero_hyperSpec$ID_Cell[i],".txt"),
        row.names=F,col.names=F,quote=F,sep = "\t")
}

#################
##normalization##
#################
data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowMeans (data_baseline_zero_hyperSpec)
#data_baseline_zero_scale_hyperSpec <- data_baseline_zero_hyperSpec / rowSums (data_baseline_zero_hyperSpec)
write.csv(data_baseline_zero_scale_hyperSpec, "Cells_bg_baseline_zero_scale.csv",
          quote=F, row.names=F)

#output txts
Cells_bg_baseline_zero_scale <- "Cells_bg_baseline_zero_scale/"
dir.create(Cells_bg_baseline_zero_scale)
for (i in (1:nrow(data_baseline_zero_scale_hyperSpec))){
  Cells<-data.frame(shift=shift,intensity=t(data_baseline_zero_scale_hyperSpec[i]$spc))
  write.table(Cells,paste0(Cells_bg_baseline_zero_scale,data_baseline_zero_scale_hyperSpec$ID_Cell[i],".txt"),
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
  Cells<-data.frame(shift=shift[1],intensity=t(cluster_means_Group[i]$spc))
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
