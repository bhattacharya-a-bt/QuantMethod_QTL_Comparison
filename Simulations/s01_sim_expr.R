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
gene_num <- as.numeric(args[1])
pass <- as.numeric(args[2])
geno_pass <- as.numeric(args[3])

# gene_num=1
# pass=1
# geno_pass=1

####################################################################################
# load dependencies silently
####################################################################################

load_silent <- function(packages) {
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE)))
  }))
}

# Usage
load_silent(c("MASS",
             "corpcor", 
             "mvnfast",
             "data.table",
             "dplyr",
             "VariantAnnotation",
             "rtracklayer",
             "Biostrings",
             "Matrix"))

####################################################################################
# helper functions
####################################################################################

# function to extract allele count from vcf genotype entry
vcf2Geno <- function(vcf){
    apply(vcf,c(1,2),FUN=function(x){
        as.numeric(substr(x,1,1))+as.numeric(substr(x,3,3))
        })
}

simulate_count <- function(gene_num) {

  set.seed(gene_num)
  print(gene_num)

  query <- anno_gene$query[gene_num]
  start <- anno_gene$start[gene_num]
  end <- anno_gene$end[gene_num]
  chr <- anno_gene$chr[gene_num]
  gene <- anno_gene$GeneID[gene_num]

  # # use query to subset the VCF with tabix
  # if(!file.exists(paste0(tmp_dir,"/subset_", chr,"_",start,"_",end, ".vcf"))){
  #   system(paste("tabix -h ", paste0(vcf_dir,"/genos_1kg_eur_500_snps_maf_0.01_chr",chr,".vcf.gz"), query, ">", 
  #              paste0(tmp_dir,"/subset_", chr,"_",start,"_",end, ".vcf")))
  # }
  
  # use query to subset the VCF with tabix
  subset_vcf_path <- paste0(tmp_dir, "/subset_", chr, "_", start, "_", end, ".vcf")

  if (!file.exists(paste0(subset_vcf_path,".gz"))) {
    system(paste("tabix -h", paste0(vcf_dir, "/genos_1kg_eur_500_snps_maf_0.01_chr", chr, ".vcf.gz"), query, ">", subset_vcf_path))
    
    # gzip and tabix the resulting VCF file
    system(paste("bgzip", subset_vcf_path))
    system(paste0("tabix -p vcf ",subset_vcf_path,".gz"))
  }

  # read in genotype data
  vcf_file <- paste0(tmp_dir,"/subset_", chr,"_",start,"_",end, ".vcf.gz")
  vcf <- readVcf(vcf_file,"hg38")
  genotypes <- geno(vcf)$GT

  # extract SNP information
  vcf_gr <- rowRanges(vcf) # Accesses the rowRanges of the VCF

  # construct SNP ID as "CHR_POS_REF_ALT"
  alt_alleles <- as.character(unlist(mcols(vcf_gr)$ALT))
  snp_ids <- paste0(seqnames(vcf_gr), "_", 
                  start(vcf_gr), "_", 
                  mcols(vcf_gr)$REF, "_", 
                  alt_alleles)


  #snp_ids <- names(vcf)

  n_snp <- nrow(genotypes)

  M <- sum(anno$GeneID == gene)  # number of isoforms
  K <- n_causal                  # number of isoQTLs per isoform
  S <- ceiling(prop_shared * K)  # number of shared isoQTLs
  D <- K - S                     # number of distinct isoQTLs
  W <- S + D * M                 # total number of causal isoQTLs
  B <- matrix(0, nrow=W, ncol=M)

  # identify causal SNPs
  causal_ind <- sample(1:n_snp, W)
  causal_shared <- causal_ind[1:S]
  X <- genotypes[causal_ind, ]
  X <- vcf2Geno(X)
  X <- t(X)
  X <- scale(X, center=TRUE, scale=TRUE)  # center and scale columns of X (??)
  N <- nrow(X)  # Nnmber of samples

  # construct cov/cor matrices for betas
  cov_B <- diag(1, M)
  cov_B[lower.tri(cov_B)] <- runif(M * (M-1) / 2, cov_B_min, cov_B_max)
  cov_B[upper.tri(cov_B)] <- t(cov_B)[upper.tri(cov_B)]
  cov_B <- cov2cor(cov_B)
  if (!is.positive.definite(cov_B)) cov_B <- make.positive.definite(cov_B)

  # shared effects (marginal var 1), then scale to per-SNP contribution
  if (S > 0) {
      h2_per_snp <- if(length(h2_g)==1) h2_g / K else (h2_g / K)  # handle scalar/vector
      shared_beta_mat <- mvrnorm(n = S, mu = rep(0, M), Sigma = cov_B)
      shared_beta_mat <- sweep(shared_beta_mat, 2, sqrt(h2_per_snp), `*`)  # column-scale
      B[1:S, ] <- shared_beta_mat
  }

  # step 2: unique effects (if any)
  if (D > 0) {
      unique_beta_vec <- rnorm(D * M, 0, sqrt(h2_g / K))
      for (m in 1:M) {
          idx_start <- S + (m - 1) * D + 1
          idx_end   <- S + m * D
          B[idx_start:idx_end, m] <- unique_beta_vec[((m - 1) * D + 1):(m * D)]
      }
  }

  # step 2: calculate the actual variance of X %*% B
  genetic_component <- X %*% B
  genetic_var <- diag(cov(genetic_component))

  # step 3: rescale B so that each isoform's genetic variance aligns with h2_g
  scaling_factors <- sqrt(h2_g / genetic_var)
  B <- B * scaling_factors

  genetic_var <- diag(cov(X %*% B))

  causal_ids <- snp_ids[causal_ind]
  #causal_ids <- paste0(causal_ind,"_",causal_ids)
  rownames(B) <- causal_ids

  # define the target variance for E term
  target_variance <- 1 - h2_g

  # generate E term covariance matrix
  cov_E <- diag(1, M)
  cov_E[upper.tri(cov_E)] <- runif(M * (M-1) / 2, cov_E_min, cov_E_max)
  cov_E[lower.tri(cov_E)] <- t(cov_E)[lower.tri(cov_E)]
  cov_E <- cov_E * target_variance / mean(diag(cov_E))
  if (!is.positive.definite(cov_E)) {
    cov_E <- make.positive.definite(cov_E)
  }
  cov_E_orig <- cov_E

  # generate expression levels (Y) for each sample
  E <- mvrnorm(n=N, mu=rep(0, M), Sigma=cov_E)
  Y <- matrix(0, nrow=N, ncol=M)

  Y <- X%*%B +E

  # calculate total variance for each isoform in Y
  total_var <- diag(cov(Y))
  
  # calculate heritability for each isoform
  heritability <- genetic_var / total_var

  transcript_ids <- anno$Transcript[anno$GeneID == gene]

  # compute gene-level expression as the sum of isoforms
  gene_expr <- rowSums(Y)

  # gene-level components
  Sigma_G <- cov(X %*% B)     # genetic covariance across isoforms
  Sigma_E <- cov(E)           # residual covariance across isoforms
  Sigma_Y <- cov(Y)           # total covariance across isoforms

  # correct multivariate gene-level h2
  one_vec <- rep(1, M)
  gene_h2 <- as.numeric( (t(one_vec) %*% Sigma_G %*% one_vec) /
                         (t(one_vec) %*% Sigma_Y %*% one_vec) )
  return(list(
    Y = Y,
    transcript_ids = transcript_ids,
    B = B,
    isoform_h2 = heritability,
    gene_h2 = gene_h2
  ))
}

####################################################################################
# generate beta matrices for isoQTLs of each gene
####################################################################################

set.seed(gene_num)

param_space <- fread(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/gene_level_parameters.txt"))
window = 250000
vcf_dir = paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",geno_pass,"/files_for_analysis/1KG_vcf")
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
#anno_gene$query <- paste(anno_gene$chr,anno_gene$window_start,anno_gene$window_end,sep=":")
anno_gene$query <- paste0(anno_gene$chr,":",anno_gene$window_start,"-",anno_gene$window_end)

gene <- anno_gene$GeneID[gene_num]
param_space <- param_space[param_space$GeneID==gene,]
n_causal = param_space$n_causal
prop_shared <- param_space$prop_shared
h2_g <- param_space$h2_g
cov_E_min = param_space$cov_E_min
cov_E_max = param_space$cov_E_max
cov_B_min = param_space$cov_B_min
cov_B_max = param_space$cov_B_max
tmp_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",geno_pass,"/files_for_analysis/tmp")

expr_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/expression")

# establish expr directory 
if(!dir.exists(expr_dir)){
  dir.create(expr_dir,recursive=T)
}

# establish tmp directory for vcf files
if(!dir.exists(tmp_dir)){
  dir.create(tmp_dir,recursive=T)
}

Y <- simulate_count(gene_num)

save(Y,file=paste0(expr_dir,"/expr_",gene,".RData"))

