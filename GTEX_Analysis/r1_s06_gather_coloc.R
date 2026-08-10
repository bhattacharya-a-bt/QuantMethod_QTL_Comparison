#!/usr/bin/env Rscript

##################################################################################
# USER CONFIGURATION - edit these paths, or set the matching environment variables
##################################################################################
# r_lib_path     extra R library path; "" to use the default .libPaths()
# coloc_out_dir  colocalization results directory written by r1_s05_coloc.R;
#                per-gene files live in '<tissue>/' and the combined hits file
#                is written at the top level
r_lib_path    <- Sys.getenv("R_LIB_PATH",    "")
coloc_out_dir <- Sys.getenv("COLOC_OUT_DIR", "/path/to/GTEx_v8/requants/coloc_results")
##################################################################################

##################################################################################
# change library to local
##################################################################################
if (nchar(r_lib_path) > 0) {
  myPaths <- .libPaths()
  myPaths <- c(r_lib_path,myPaths)
  .libPaths(myPaths)
}

####################################################################################
# load dependencies
####################################################################################
library(data.table)
library(dplyr)
library(stringr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
annot <- as.character(args[1])
quant <- as.character(args[2])
tissue <- as.character(args[3])

####################################################################################
# helper functions
####################################################################################

pheno_name_vec <- c("BreastCancer",
  "ProstateCancer",
  "T2Diabetes",
  "BMI",
  "Height",
  "Schizophrenia",
  "BipolarDisorder")

setwd(file.path(coloc_out_dir,tissue))

files <- list.files()

for(pheno_name in pheno_name_vec){

print(pheno_name)

foc <- files[grep(pheno_name,files)]
foc <- foc[grep(annot,foc)]
foc <- foc[grep(quant,foc)]

if(length(foc)>0){

dat <- {}
for(i in 1:length(foc)){
  #print(i)
  tmp <- data.frame(fread(foc[i]))
  tmp <- subset(tmp, CLPP >= 0.01 & P <= 5e-7)
  if(nrow(tmp)>0){
    dat <- rbind(dat,tmp)
  }
}

if (nrow(dat) > 0) {
  cat('Writing out!\n')
  
  dat$quant <- quant
  dat$annot <- annot
  dat$tissue <- tissue
  dat$pheno <- pheno_name

  # create unique output filename for this tissue-cancer combination
  output_filename <- paste0('colocalization_results_', tissue, '_', annot,"_",quant,"_",pheno_name, '.tsv')
  
  fwrite(dat,
         file.path(coloc_out_dir, output_filename),
         row.names = FALSE,
         quote = FALSE,
         sep = '\t')
  
  cat(paste('Results written to:', output_filename, '\n'))
} else {
  cat('No hits\n')
}

cat(paste('Completed processing tissue', tissue, 'cancer', pheno_name, '\n'))

}else{
  print("No matching colocalization files for these annot-quant-tissue-pheno combination.")
} 

}

