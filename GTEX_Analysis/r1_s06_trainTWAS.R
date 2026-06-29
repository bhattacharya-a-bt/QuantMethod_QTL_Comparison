#!/usr/bin/env Rscript

##################################################################################
# change library to local
##################################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

####################################################################################
# load dependencies
####################################################################################
library(data.table)
library(isotwas)
library(bigsnpr)
library(dplyr)
library(stringr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
annot <- as.character(args[1])
quant <- as.character(args[2])
tissue <- as.character(args[3])
index <- as.integer(args[4])

# check if index conversion was successful
if (is.na(index)) {
  stop("The fourth argument (index) must be a valid numeric value.", call. = FALSE)
}

PSPACE <- data.frame(fread("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt",
  header=F))

cat("Processing tissue:", tissue, "\n")
cat("Processing annot:", annot, "\n")
cat("Processing quant:", quant, "\n")
cat("Index:", index, "\n")

bed_file = paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/",tissue,"/",annot,"_",quant,".v8.normalized_expression.bed.gz")
bed <- data.frame(fread(bed_file))

temp_dir = file.path('/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/tmp',tissue)
dir.create(temp_dir,recursive = T)
VCFFILE="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated.vcf.gz"
COVARFILE=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/",tissue,"/",tissue,"_formatted_covariates.txt")

cat('Covariate residualizing')
outfile = paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/",tissue,"/",annot,"_",quant,".v8.normalized_expression_corrected.bed")
if (!file.exists(outfile)){
  system(paste('QTLtools correct --bed ',bed_file,
               '--cov',COVARFILE,'--normal --out',outfile))}

bed <- data.frame(fread(outfile))

gtex_snpfile = '/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated'
snps = snp_attach(snp_readBed2(paste0(gtex_snpfile,'.bed'),
                               backingfile = tempfile()))
dir.create(temp_dir,recursive = T)
snps_matrix = as.matrix(snps$genotypes[])
colnames(snps_matrix) = snps$map$marker.ID
rownames(snps_matrix) = snps$fam$sample.ID

result_df <- data.frame(tissue=tissue,
  version=annot, method=quant)

# function to process TWAS models for a single gene (one row)
process_single_gene <- function(row_index, gene_data, snps_matrix, snps) {
  
  out_dir = file.path('/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/TWAS_weights',
                      paste(tissue,annot,quant,sep="_"))
  dir.create(out_dir, recursive = T, showWarnings = FALSE)
  
  cat('Extracting information for row', row_index, '...\n')
  colnames(gene_data)[1] <- "#chr"
  gene = gene_data$id[row_index]
  chr = gene_data$`#chr`[row_index]
  chr <- as.integer(substr(chr,4,nchar(chr)))
  start = as.integer(gene_data$start[row_index])
  end = as.integer(gene_data$end[row_index])
  
  if (!file.exists(file.path(out_dir, paste0(gene, '_TWAS.RDS')))){
    
    samps_in_bed <- colnames(gene_data)[7:ncol(gene_data)]
    samps_in_geno <- rownames(snps_matrix)
    samps_keep <- intersect(samps_in_bed, samps_in_geno)

    gene_data <- gene_data[, c(colnames(gene_data)[1:6], samps_keep)]

    snps_matrix_this = snps_matrix[colnames(gene_data)[7:ncol(gene_data)],
                                   snps$map$marker.ID[snps$map$chromosome == chr &
                                                        snps$map$physical.pos < end + 1e6 &
                                                        snps$map$physical.pos > start - 1e6]]

    stopifnot(identical(colnames(gene_data)[7:ncol(gene_data)],
                    rownames(snps_matrix_this)))
    
    gene_exp_vector = as.numeric(as.data.frame(gene_data)[row_index, 7:ncol(gene_data)])
    
    cat(paste0('Training models for ', gene, ' at Chr', chr, ' ', start, ':', end, '\n'))
    
    enet = univariate_elasticnet(X = snps_matrix_this,
                                 Y = as.matrix(gene_exp_vector),
                                 Omega = NA,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = NULL,
                                 seed = NULL)
    
    blup = univariate_blup(X = snps_matrix_this,
                           Y = as.matrix(gene_exp_vector),
                           Omega = NA,
                           scale = F,
                           nfolds = 5,
                           verbose = F,
                           par = F,
                           n.cores = NULL,
                           tx_names = NULL,
                           seed = NULL)
    
    susie = univariate_susie(X = snps_matrix_this,
                             Y = as.matrix(gene_exp_vector),
                             Omega = NA,
                             scale = F,
                             alpha = 0.5,
                             nfolds = 5,
                             verbose = F,
                             par = F,
                             n.cores = NULL,
                             tx_names = NULL,
                             seed = NULL)
    
    cat('Finding optimal model\n')
    r2 = unlist(c(enet[[1]]$R2,
                  blup[[1]]$R2,
                  susie[[1]]$R2))
    
    if (max(r2) > .01){
      r2.min = which(r2 == max(r2))
      if (r2.min == 1){
        last.mod = enet
        if (nrow(enet[[1]]$Model) == 0){
          r2.min = which(r2 == median(r2))
        }
      }
      
      if (r2.min == 2){
        last.mod = blup
      }
      
      if (r2.min == 3){
        last.mod = susie
      }
      
      gene.model = last.mod[[1]]$Model
      gene.model$Feature = gene
      gene.model$R2 = last.mod[[1]]$R2
      gene.model$marker.ID = gene.model$SNP
      gene.model = merge(gene.model, snps$map, by='marker.ID')
      colnames(gene.model) = c('marker.ID','SNP',
                               'Weight','Feature','R2','Chromosome',
                               'genetic.dist','Position',
                               'ALT','REF')
      gene.model$Build = 'hg38'
      gene.model = gene.model[,c('Feature','SNP','Chromosome','Position','Build',
                                 'ALT','REF','Weight','R2')]
      
      saveRDS(gene.model, file.path(out_dir, paste0(gene, '_TWAS.RDS')))
    } else {
      saveRDS(paste0('R2 < 0.01'), file.path(out_dir, paste0(gene, '_TWAS.RDS')))
    }
  }
  
  return(gene) # return gene name for tracking
}

# function to process all rows for a given dataset using lapply
process_twas_models <- function(gene_data, snps_matrix, snps) {
  cat("Processing", nrow(gene_data), "genes...\n")
  
  # use lapply to process each row with error handling
  results <- lapply(1:nrow(gene_data), function(i) {
    tryCatch({
      process_single_gene(i, gene_data, snps_matrix, snps)
    }, error = function(e) {
      cat("Error processing gene at row", i, "in:", e$message, "\n")
      return(paste0("ERROR_ROW_", i))
    })
  })
  
  # count successful vs failed processes
  errors <- sum(grepl("^ERROR_ROW_", results))
  successes <- length(results) - errors
  
  cat("Completed processing:", successes, "successful,", errors, "errors\n")
  return(results)
}

# bin gene list in bed files
set.seed(123)
n_bins=10
gene_count <- 1:nrow(bed)
bins <- cut(
  gene_count,
  breaks = n_bins,
  labels = FALSE
)
gene_data = bed[bins==index,]
gene_data$start <- as.numeric(gene_data$start)
gene_data$end <- as.numeric(gene_data$end)

# process TWAS models
process_twas_models(gene_data, snps_matrix, snps)

