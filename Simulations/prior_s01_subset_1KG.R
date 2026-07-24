#!/usr/bin/env Rscript

###################################################################
# edit these paths and parameters before running
###################################################################

## ---- static parameters ----
pass <- 2 # run apss
n <- 500 # number of subjects to extract from 1KG for simulation
seed <- 12345
pop <- "EUR" # 1KG population to simulation from

## ---- input reference file ----
geno_file <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/all_hg38.psam"

## ---- project / working directory ----
# project directory base; the pass number is appended automatically.
# kgp3.RData is read from <proj_dir>/files_for_analysis/, from kgp R package 
proj_base_dir <- "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp"

## ---- output directory ----
# output directory base; the pass number is appended automatically
out_base_dir <- "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp"

####################################################################################
# derived paths (do not edit below unless changing structure)
####################################################################################
proj_dir <- paste0(proj_base_dir, "/pass", pass)
directory <- paste0(out_base_dir, "/pass", pass, "/files_for_analysis")

####################################################################################
# setup output directory
####################################################################################

# check if the directory doesn't exist
if (!file.exists(directory)) {
  # Create the directory along with any intermediate directories
  dir.create(directory, recursive = TRUE)
  cat("Directory created:", directory, "\n")
} else {
  cat("Directory already exists:", directory, "\n")
}

setwd(proj_dir)

####################################################################################
# select EUR samples from 1KG
####################################################################################

dat <- read.csv(geno_file, sep = "\t")
eur <- dat[dat$SuperPop == pop, ]
table(eur$Population)
dim(eur) # 633x6

# read in 1KG metadata from "kgp" R package
# cannot install this package on HPC due to R version
load("files_for_analysis/kgp3.RData")
dim(kgp3) # 2504 samples sequenced for the phase 3 release

eur <- eur[eur$X.IID %in% kgp3$id, ]
dim(eur) # 503x6

# randomly select n
set.seed(seed)
sel <- sample(1:nrow(eur), n)
samp <- data.frame(ID = eur[sel, 1])

write.table(samp, file = "files_for_analysis/1kg_eur_500_sample_ids",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
