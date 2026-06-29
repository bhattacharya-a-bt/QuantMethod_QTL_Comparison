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
library(rtracklayer)
library(GenomicRanges)
library(bigsnpr)
library(stringr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
index <- as.integer(args[1])
gwas <- as.character(args[2])
pheno <- as.character(args[3])
liftoverTo38 <- as.logical(args[4])
nbins <- as.integer(args[5])
bin <- as.integer(args[6])

run_dat <- readRDS("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/twas_results/r1_aggregated_TWAS_passR2.RDS")

run_dat$annot[run_dat$annot=="GENCODE_V27"] <- "GENCODE_v27"
run_dat$annot[run_dat$annot=="GENCODE_V38"] <- "GENCODE_v38"
run_dat$annot[run_dat$annot=="GENCODE_V45"] <- "GENCODE_v45"

run_dat <- run_dat[index,]
gene_list <- run_dat$pass_genes[1][[1]][[1]]
gwas_name <- pheno

# calculate chunk for distributing jobs into 20 arrays
total_files <- length(gene_list)
chunk <- ceiling(total_files / nbins)
index_start <- (bin - 1) * chunk + 1
index_stop <- min(index_start + chunk - 1, total_files)

# GWAS file with their names for output
gwas_file <- file.path('/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS',
                   c(gwas))

# chain file for liftover
chain_file <- '/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain'

# function to perform liftover using R liftOver package
liftover_gwas <- function(gwas_data, chain_file) {
  require(rtracklayer)
  require(GenomicRanges)
  require(data.table)
  
  # load the chain file
  chain <- import.chain(chain_file)
  
  # create GRanges object
  gwas_gr <- GRanges(
    seqnames = paste0("chr", gwas_data$CHR),
    ranges = IRanges(start = gwas_data$POS, end = gwas_data$POS),
    strand = "*"
  )
  
  # add metadata to track original rows
  mcols(gwas_gr) <- data.frame(
    orig_idx = 1:nrow(gwas_data),
    SNP = gwas_data$SNP,
    CHR = gwas_data$CHR,
    POS = gwas_data$POS
  )
  
  # liftover
  lifted_gr <- liftOver(gwas_gr, chain)
  
  # extract successfully lifted variants
  lifted_gr_unlisted <- unlist(lifted_gr)
  
  if (length(lifted_gr_unlisted) > 0) {
    # get  original indices of lifted variants
    orig_indices <- mcols(lifted_gr_unlisted)$orig_idx
    
    gwas_lifted <- gwas_data[orig_indices, ]
    
    # update coordinates with lifted positions
    gwas_lifted$CHR <- as.numeric(gsub("chr", "", as.character(seqnames(lifted_gr_unlisted))))
    gwas_lifted$POS <- start(lifted_gr_unlisted)
    
    cat("    Liftover successful:", length(lifted_gr_unlisted), "out of", nrow(gwas_data), "variants lifted\n")
    
    return(gwas_lifted)
  } else {
    warning("Liftover failed for all variants, using original coordinates")
    return(gwas_data)
  }
}

# GTEx SNP reference file
gtex_snpfile = '/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated'

cat("Processing GWAS:", gwas_name, "\n")
cat("Processing array index:", bin, "covering files", index_start, "to", index_stop, "\n")

# load and liftover gwas data
gwas_data <- data.table::fread(gwas_file)

# check if columns match expected format and rename if necessary
if(gwas_name %in% c("BreastCancer","ProstateCancer")){
  if ("SNP" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
  # rename columns to match expected format
  colnames(gwas_data)[colnames(gwas_data) == "BP"] <- "POS"
  colnames(gwas_data)[colnames(gwas_data) == "A1"] <- "ALT"
  colnames(gwas_data)[colnames(gwas_data) == "A2"] <- "REF"
  } else {
    warning("GWAS file format may not match expected columns")
  }
  }else{

if ("SNP" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
    # Rename columns to match expected format
    colnames(gwas_data)[colnames(gwas_data) == "BP"] <- "POS"
    colnames(gwas_data)[colnames(gwas_data) == "A1"] <- "REF"
    colnames(gwas_data)[colnames(gwas_data) == "A2"] <- "ALT"
  } else {
    warning("GWAS file format may not match expected columns")
  }
}

# remove prefix chr if present
if(length(grep("^chr",gwas_data$CHR[1]))==1){
    gwas_data[, CHR := sub("^chr", "", CHR)]
    gwas_data$CHR <- as.numeric(gwas_data$CHR)
}

# liftover on full GWAS data if needed
if(liftoverTo38==T){
  gwas_data_lifted <- liftover_gwas(gwas_data,chain_file)
}else{
  gwas_data_lifted <- gwas_data
}
   
# calculate Z-score if not present
if (!"Z" %in% colnames(gwas_data_lifted)) {
  gwas_data_lifted$Z <- gwas_data_lifted$BETA / gwas_data_lifted$SE
}

gwas_data_lifted$CHRPOS <- paste(gwas_data_lifted$CHR, gwas_data_lifted$POS, sep = ':')

results <- data.frame()

annot=run_dat$annot
quant=run_dat$quant
tissue=run_dat$tissue

# loop through each TWAS model file
for (i in index_start:index_stop) {
  
  cat("Processing TWAS model", index, "of", total_files, "\n")
  
  gene <- gene_list[i]
  twas_file <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/TWAS_weights/",tissue,"_",annot,"_",quant,"/",gene,"_TWAS.RDS")  
  # Check if file exists
  if (!file.exists(twas_file)) {
    cat("  TWAS file not found:", twas_file, "\n")
    next
  }
  
  # read TWAS model
  twas_model <- readRDS(twas_file)
  
  if (nrow(twas_model) == 0) {
    cat("  Empty TWAS model\n")
    next
  }
  
  chr <- twas_model$Chromosome[1]
  min_pos <- min(twas_model$Position)
  max_pos <- max(twas_model$Position)
  
  cat("  Gene:", gene, "Tissue:", tissue, "Version:", annot, "Method:", quant, "\n")
  
  twas_model$CHRPOS <- paste(twas_model$Chromosome, twas_model$Position, sep = ':')
  
  # create temporary folder
  tempfolder <- file.path('/rsrch5/home/epi/sthead/',
                          'tempTWAS',
                          paste0('gene', i, '_', gene, '_', gwas_name, "_",tissue,"_",annot,"_",quant))
  dir.create(tempfolder, recursive = TRUE, showWarnings = FALSE)
  
  # use plink to extract SNPs for this region
  system(paste('plink2 --bfile', gtex_snpfile,
               '--chr', chr,
               '--from-bp', min_pos,
               '--to-bp', max_pos,
               '--make-bed',
               '--out', file.path(tempfolder, paste0(gene, '_', tissue))))
  
  bed_file_path <- file.path(tempfolder, paste0(gene, '_', tissue, '.bed'))
  
  if (!file.exists(bed_file_path)) {
    cat("  Failed to create bed file\n")
    unlink(tempfolder, recursive = TRUE)
    next
  }
  
  # read genotype data and compute LD
  tryCatch({
    snps <- snp_attach(snp_readBed2(bed_file_path,
                                    backingfile = file.path(tempfolder, paste0(gene, '_', tissue, '.bk'))))
    
    snps$map$CHRPOS <- paste(snps$map$chromosome, snps$map$physical.pos, sep = ':')
    
    # compute LD matrix
    ld <- snp_cor(snps$genotypes)
    colnames(ld) <- rownames(ld) <- snps$map$CHRPOS
    
    # find intersection of GWAS and TWAS model SNPs
    int_chrpos <- intersect(gwas_data_lifted$CHRPOS, twas_model$CHRPOS)
    int_chrpos <- intersect(int_chrpos, snps$map$CHRPOS)
    
    if (length(int_chrpos) == 0) {
      cat("  No overlapping SNPs between GWAS, TWAS model, and genotype data\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    # sub data to overlapping SNPs
    gwas_subset <- gwas_data_lifted[gwas_data_lifted$CHRPOS %in% int_chrpos, ]
    twas_subset <- twas_model[twas_model$CHRPOS %in% int_chrpos, ]
    snps_subset_idx <- which(snps$map$CHRPOS %in% int_chrpos)
    ld_subset <- ld[int_chrpos, int_chrpos]
    
    gwas_subset <- gwas_subset[match(int_chrpos, gwas_subset$CHRPOS), ]
    twas_subset <- twas_subset[match(int_chrpos, twas_subset$CHRPOS), ]
    snps_map_subset <- snps$map[snps_subset_idx, ][match(int_chrpos, snps$map[snps_subset_idx, ]$CHRPOS), ]
    
    # check and flip alleles if necessary
    twas_weights <- twas_subset$Weight
    
    # compare ALT alleles between GWAS and TWAS model
    for (i in 1:length(int_chrpos)) {
      if (!is.na(gwas_subset$ALT[i]) && !is.na(snps_map_subset$allele1[i])) {
        if (gwas_subset$ALT[i] != snps_map_subset$allele1[i]) {
          # flip the weight if alleles don't match
          twas_weights[i] <- -twas_weights[i]
        }
      }
    }
    
    # compute TWAS z-statistic
    # formula: (Weight * GWAS_Z) / sqrt(Weight^T * LD * Weight)
    gwas_z <- gwas_subset$Z
    
    numerator <- sum(twas_weights * gwas_z)
    denominator <- sqrt(as.numeric(t(twas_weights) %*% ld_subset %*% twas_weights))
    
    if (denominator == 0 || is.na(denominator)) {
      cat("  Invalid denominator for TWAS Z-statistic\n")
      twas_z <- NA
    } else {
      twas_z <- numerator / denominator
    }
    
    # find top GWAS SNP within 1e6 bases
    extended_region <- gwas_data_lifted[gwas_data_lifted$CHR == chr & 
                                          gwas_data_lifted$POS >= (min_pos - 1e6) &
                                          gwas_data_lifted$POS <= (max_pos + 1e6), ]
    
    if (nrow(extended_region) > 0) {
      top_snp_idx <- which.min(extended_region$P)
      top_snp_pval <- extended_region$P[top_snp_idx]
      top_snp_id <- extended_region$SNP[top_snp_idx]
    } else {
      top_snp_pval <- NA
      top_snp_id <- NA
    }
    
    # add result to data frame
    result_row <- data.frame(
      Cancer = gwas_name,
      Gene = gene,
      Tissue = tissue,
      Version = annot,
      Method = quant,
      TWAS_Z = twas_z,
      Top_GWAS_SNP = top_snp_id,
      Top_GWAS_P = top_snp_pval,
      N_SNPs = length(int_chrpos),
      stringsAsFactors = FALSE
    )
    
    
    # write out individual result
    dir_out <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/twas_results/",tissue,"/",annot)
    dir.create(dir_out, recursive = TRUE)
    
    data.table::fwrite(result_row,
                               file.path(dir_out, paste0(gene, '_twas_',gwas_name,'_',annot,'_',quant,'.txt.gz')),
                               sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    cat("  TWAS Z-score:", twas_z, "\n")
    cat('Written out to',
        file.path(dir_out, paste0(gene, '_eQTL_coloc_',gwas_name,'_',annot,'_',quant,'.txt.gz')))
    
    
  }, error = function(e) {
    cat("  Error processing:", e$message, "\n")
  })
  
  # clean up temporary files
  unlink(tempfolder, recursive = TRUE)
}

cat("Completed processing", nrow(results), "TWAS models for", gwas_name, "\n")

