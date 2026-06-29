#!/usr/bin/env Rscript

####################################################################################
# change library to local
####################################################################################

myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3", myPaths)
.libPaths(myPaths)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
sample <- as.numeric(args[1])
pass <- as.numeric(args[2])
param_row_reads <- as.numeric(args[3])

####################################################################################
# load dependencies
####################################################################################
library(MASS)
library(corpcor)
library(mvnfast)
library(data.table)
library(Rsubread)
library(dplyr)
library(VariantAnnotation)
library(rtracklayer)
library(Biostrings)
library(stringr)

####################################################################################
# set dummy data
####################################################################################
#gene_num = 1
window = 250000

####################################################################################
# read in annotation data
####################################################################################

anno <- read.table(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/anno_selected_genes.txt"),header=T)

# summarize annotation table by gene
anno_gene <- anno %>%
  group_by(GeneID) %>%
  summarise(start = unique(start), end = unique(end), seqnames = unique(seqnames), strand = unique(strand))

# format the coordinates for tabix
anno_gene$window_start <- anno_gene$start - window
anno_gene$window_start[anno_gene$window_start<1] <- 1
anno_gene$window_end <- anno_gene$end + window
anno_gene$chr <- as.numeric(substr(anno_gene$seqnames,4,nchar(anno_gene$seqnames)))
anno_gene$query <- paste(anno_gene$chr,anno_gene$window_start,anno_gene$window_end,sep=":")

gene_list <- unique(anno_gene$GeneID)

####################################################################################
# simulate reads
####################################################################################

param_space_reads <- fread(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/parameter_space_reads.txt"))
library_size <- param_space_reads$library_size[param_row_reads]
read_length <- param_space_reads$read_length[param_row_reads]
paired_end <- param_space_reads$paired_end[param_row_reads]
paired_end <- as.logical(paired_end)

reads_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/reads")

if(!dir.exists(reads_dir)){
  dir.create(reads_dir)
}

fa.file <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.transcripts.fa")
transcripts <- scanFasta(fa.file)
tx_id <- str_extract(transcripts$TranscriptID, "^[^|]+")


if(file.exists(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/expression/aggregatedY.RData"))){
   load(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/expression/aggregatedY.RData"))
}else{

  Y_files <- list.files(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/expression"),full=T)

  Y_list <- lapply(Y_files, function(file) {
    load(file)
    return(Y)
  })

  Y_full <- do.call(cbind, lapply(Y_list, `[[`, 1))  # combine all Y[[1]] into a matrix
  tx_full <- unlist(lapply(Y_list, `[[`, 2))        # combine all Y[[2]] into a vector
  save(Y_full,tx_full,file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/expression/aggregatedY.RData"))

}

Y <- data.frame(expr=Y_full[sample,],
                tx=str_extract(tx_full, "^[^|]+"))

Y$expr <- exp(Y$expr)
id <- rownames(Y_full)[sample]

TPMs <- data.frame(tx=tx_id)
TPMs$id <- 1:nrow(TPMs)

dat <- merge(x=TPMs,y=Y,by="tx",all.x=T,all.y=F)
dat$expr[is.na(dat$expr)] <- 0
dat <- dat[order(dat$id),]
sum(dat$tx==tx_id)

setwd(reads_dir)

simReads(transcript.file = paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.transcripts.fa"),
         expression.levels=dat$expr,
         output.prefix=paste0('sim_',id,"_param_row_reads_",param_row_reads),
         read.length = read_length,         
         library.size = library_size,     
         paired.end = paired_end)    


