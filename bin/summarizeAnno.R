#!/usr/bin/env Rscript

## summarize gene annotations and contig annotations
# wget http://rest.kegg.jp/list/ko -O ko_def.txt

setwd('./')

# Load / install packages
if (!require("ape")) {
  install.packages("ape", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("ape")
}

library(optparse)

option_list = list(
    make_option(c("-c", "--cell"), type="character", default="RG1", help="cell prefix", metavar="character"),
    make_option(c("-K", "--KO"), type="character", default="ko_KO.txt", help="ko-KO definition file name [default= %default]", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="RG1_annotations.txt", help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ko_KO <- read.delim(opt$KO, header=F, stringsAsFactors=F, colClasses="character")

gff_file <- paste0("../results/prokka/",opt$cell,"/",opt$cell,".gff")

raw_gff <- readLines(gff_file)
fasta_tag_line <- grep("^##FASTA", raw_gff)
new_gff <- raw_gff[1:fasta_tag_line-1]

gff_path <- dirname(gff_file)
gff_name <- basename(gff_file)
new_name <- paste0(gsub(pattern = "\\.gff$", "", gff_name),"_fix.gff3")
new_file <- paste(gff_path,new_name,sep="/")

ok <- file.create(new_file)
writeLines(new_gff, new_file)

gene_model <- read.gff(new_file)
genes <- gene_model[gene_model$type=="gene",]
gids <- sapply(strsplit(genes$attributes, "locus_tag="), function(item) {item[2]})

g2c <- data.frame(ctg=genes$seqid, gid=gids, stringsAsFactors=F)

blob_file <- paste0("../results/blob/",opt$cell,"/",opt$cell,".blobDB.bestsum.table.txt")
blob_tbl <- read.delim(blob_file, header=F, comment.char="#", stringsAsFactors=F, row.names=NULL)
ctg_ids <- sapply(strsplit(blob_tbl$V1, "_length_"), function(item) {item[1]})
blob_tbl$ctg <- ctg_ids
names(blob_tbl)[names(blob_tbl) == "V21"] <- "Species"

g2c2s <- merge(g2c, blob_tbl[,c("ctg","Species")], by="ctg", all.x=T)

egg_file <- paste0("../results/eggnog/",opt$cell,".emapper.annotations")
egg_nog <- read.delim(egg_file, header=F, comment.char="#", stringsAsFactors=F, row.names=NULL)
names(egg_nog)[names(egg_nog) == "V1"] <- "gid"
names(egg_nog)[names(egg_nog) == "V7"] <- "KOs"

eggnog <- egg_nog[egg_nog$KOs!="",c("gid","KOs")]

data <- merge(g2c2s, eggnog, by="gid", all.x=T)

ko_KOs <- sapply(data$KOs, function(item) {
    if(!is.na(item)) {
        KOs <- strsplit(item, ",")
        kos <- c()
        for (KO in KOs) {
            ko <- ko_KO[ko_KO$V3==KO,]
            kos <- c(kos, paste0("ko",ko$V1," ",ko$V2," ||| ",KO," ",ko$V4))
        }
        paste(kos,collapse = ' +++ ')
    } else { NA }
})

data$KO_def <- ko_KOs

write.table(data, file=opt$out, append=F, quote=F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names=F, col.names=T)
