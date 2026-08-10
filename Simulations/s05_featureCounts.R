#!/usr/bin/env Rscript

####################################################################################
# USER CONFIGURATION - edit these paths, or set the matching environment variables
####################################################################################
# r_lib_path  extra R library path; "" to use the default .libPaths()
# base_dir    base project directory; the pass number is appended internally
r_lib_path <- Sys.getenv("R_LIB_PATH", "")
base_dir   <- Sys.getenv("BASE_DIR",   "/path/to/scratch/GTEx_gencode_comp")
####################################################################################

####################################################################################
# change library to local
####################################################################################
if (nchar(r_lib_path) > 0) {
  myPaths <- .libPaths()
  myPaths <- c(r_lib_path, myPaths)
  .libPaths(myPaths)
}

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
pass <- as.numeric(args[1])
geno_pass <- as.numeric(args[2])
param_row_reads <- as.numeric(args[3])
gtf <- as.character(args[4])
annot_name <- as.character(args[5])

####################################################################################
# load dependencies silently
####################################################################################
suppressPackageStartupMessages({
  library(Rsubread)
  library(rtracklayer)
  library(data.table)
  library(optparse)
})

SAMPLES_FILE=paste0(base_dir,"/pass",geno_pass,"/files_for_analysis/1kg_eur_500_sample_ids")
samples = data.table::fread(SAMPLES_FILE,header=F)$V1

# determine whether single or paired end reads
dat <- fread(paste0(base_dir,"/pass",pass,"/files_for_analysis/parameter_space_reads.txt"))
paired_end <- dat$paired_end[param_row_reads]

## annotation: GeneID Chr Start  End Strand
gtf.df = data.frame(import(gtf))
annot = gtf.df[gtf.df$type=="gene",]
annot = annot[,c("gene_id","seqnames","start","end","strand")]
names(annot) = c("GeneID","Chr","Start","End","Strand")

DIR_BAM=paste0(base_dir,"/pass",pass,"/files_for_analysis/star_alignments/param_row_reads_",param_row_reads)

bams = file.path(DIR_BAM,paste0(samples,"_Aligned.sortedByCoord.out.bam"))

if(paired_end=="T"){
  counts = featureCounts(bams,annot.ext=annot,isPairedEnd = TRUE, nthreads=10) # check documentation RE: multi-mapping and multi-feature mapping
}else{
  counts = featureCounts(bams,annot.ext=annot,isPairedEnd = FALSE, nthreads=10) # check documentation RE: multi-mapping and multi-feature mapping
}

saveRDS(list(gtf=gtf,annot = annot,
             counts = counts),
        file.path(DIR_BAM,
          paste0(annot_name,'.RDS')))



