#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(bigsnpr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript script.R <array_index>")
}

gene_info <- readRDS('/rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/eGene_colocalization.RDS')
array_index <- as.numeric(args[1])
chunk <- 25
index_start <- (array_index - 1) * chunk + 1
index_stop <- index_start + chunk - 1
index_stop = min(index_stop, length(unique(gene_info$ensembl_gene_id)))

# GWAS files with their names for output
gwas_files <- data.frame(
  path = file.path('/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS',
                   c('Endometrial.tsv.gz',
                     'Renal.tsv.gz',
                     'ProstateOverall.tsv.gz',
                     'OvarianOverall.tsv.gz',
                     'Melanoma.tsv.gz',
                     'LungOverall.tsv.gz',
                     'HeadNeck.tsv.gz',
                     'Glioma.tsv.gz',
                     'Colorectal.tsv.gz',
                     'BreastOverall.tsv.gz')),
  name = c('Endometrial', 'Renal', 'ProstateOverall', 'OvarianOverall',
           'Melanoma', 'LungOverall', 'HeadNeck', 'Glioma', 'Colorectal', 'BreastOverall'),
  stringsAsFactors = FALSE
)

# Chain file for liftover
chain_file <- '/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain'

# Function to perform liftover using R liftOver package
liftover_gwas <- function(gwas_data, chain_file) {
  require(rtracklayer)
  require(GenomicRanges)
  require(data.table)
  
  # Load the chain file
  chain <- import.chain(chain_file)
  
  # Create GRanges object from GWAS data
  gwas_gr <- GRanges(
    seqnames = paste0("chr", gwas_data$CHR),
    ranges = IRanges(start = gwas_data$POS, end = gwas_data$POS),
    strand = "*"
  )
  
  # Add metadata to track original rows
  mcols(gwas_gr) <- data.frame(
    orig_idx = 1:nrow(gwas_data),
    SNP = gwas_data$SNP,
    CHR = gwas_data$CHR,
    POS = gwas_data$POS
  )
  
  # Perform liftover
  lifted_gr <- liftOver(gwas_gr, chain)
  
  # Extract successfully lifted variants
  lifted_gr_unlisted <- unlist(lifted_gr)
  
  if (length(lifted_gr_unlisted) > 0) {
    # Get the original indices of successfully lifted variants
    orig_indices <- mcols(lifted_gr_unlisted)$orig_idx
    
    # Create lifted GWAS data
    gwas_lifted <- gwas_data[orig_indices, ]
    
    # Update coordinates with lifted positions
    gwas_lifted$CHR <- as.numeric(gsub("chr", "", as.character(seqnames(lifted_gr_unlisted))))
    gwas_lifted$POS <- start(lifted_gr_unlisted)
    
    cat("    Liftover successful:", length(lifted_gr_unlisted), "out of", nrow(gwas_data), "variants lifted\n")
    
    return(gwas_lifted)
  } else {
    warning("Liftover failed for all variants, using original coordinates")
    return(gwas_data)
  }
}

# Loop through each gene
for (index in index_start:index_stop) {
  
  gene_interest <- subset(gene_info, ensembl_gene_id == unique(gene_info$ensembl_gene_id)[index])
  
  gene <- gene_interest$ensembl_gene_id[1]
  start <- gene_interest$start_position[1]
  chr <- gene_interest$chromosome_name[1]
  end <- gene_interest$end_position[1]
  biotype <- gene_interest$gene_biotype[1]
  
  cat("Processing gene:", gene, "on chromosome", chr, "\n")
  
  # Loop through each GWAS file
  for (gwas_idx in 1:nrow(gwas_files)) {
    
    gwas_file <- gwas_files$path[gwas_idx]
    gwas_name <- gwas_files$name[gwas_idx]
    
    cat("  Processing GWAS:", gwas_name, "\n")
    
    
    # Create temp folder for this gene-GWAS combination
    tempfolder <- file.path('/rsrch5/scratch/epi/abhattacharya3',
                            'tempColoc',
                            paste0('gene', index, '_', gene, '_', gwas_name))
    dir.create(tempfolder, recursive = TRUE)
    
    ### LOAD AND LIFTOVER GWAS DATA
    gwas_data <- data.table::fread(gwas_file)
    
    # Check if columns match expected format and rename if necessary
    if ("SNP" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
      # Rename columns to match expected format
      colnames(gwas_data)[colnames(gwas_data) == "BP"] <- "POS"
      colnames(gwas_data)[colnames(gwas_data) == "A1"] <- "ALT"
      colnames(gwas_data)[colnames(gwas_data) == "A2"] <- "REF"
    } else {
      warning("GWAS file format may not match expected columns")
    }
    
    # Filter to chromosome region (pre-liftover)
    gwas_data_subset <- subset(gwas_data,
                               CHR == chr)
    
    if (nrow(gwas_data_subset) == 0) {
      cat("    No SNPs found in region for", gwas_name, "\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    # Perform liftover
    gwas_data_lifted <- liftover_gwas(gwas_data_subset, chain_file)
    
    # Calculate Z-score if not present
    if (!"Z" %in% colnames(gwas_data_lifted)) {
      gwas_data_lifted$Z <- gwas_data_lifted$BETA / gwas_data_lifted$SE
    }
    
    # Filter again after liftover to hg38 coordinates
    gwas_data_lifted <- subset(gwas_data_lifted,
                               CHR == chr &
                                 POS < end + 1e6 &
                                 POS > start - 1e6)
    
    if (nrow(gwas_data_lifted) == 0 | min(gwas_data_lifted$P) <= 1e-7) {
      cat("    No SNPs found in region after liftover for", gwas_name, "\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    gwas_data_lifted$CHRPOS <- paste(gwas_data_lifted$CHR, gwas_data_lifted$POS, sep = ':')
    
    
    ### EXTRACT LD FROM GTEx
    ecaviar <- '/rsrch5/home/epi/bhattacharya_lab/software/caviar/CAVIAR-C++/eCAVIAR'
    gtex_snpfile <- '/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF'
    
    require(bigsnpr)
    system(paste('plink2 --bfile', paste0(gtex_snpfile),
                 '--chr', chr,
                 '--from-bp', start - 1e6,
                 '--to-bp', end + 1e6,
                 '--make-bed',
                 '--out', file.path(tempfolder, paste0(gene, '_', gwas_name))))
    
    bed_file_path <- file.path(tempfolder, paste0(gene, '_', gwas_name, '.bed'))
    if (!file.exists(bed_file_path)) {
      cat("    Failed to create GTEx bed file for", gwas_name, "\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    snps_gtex <- snp_attach(snp_readBed2(bed_file_path))
    snps_gtex$map$CHRPOS <- paste(snps_gtex$map$chromosome,
                                  snps_gtex$map$physical.pos,
                                  sep = ':')
    
    ### INTERSECT GWAS AND GTEX
    int_chrpos <- intersect(gwas_data_lifted$CHRPOS, snps_gtex$map$CHRPOS)
    
    if (length(int_chrpos) == 0) {
      cat("    No overlapping SNPs between GWAS and GTEx for", gwas_name, "\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    gwas_data_lifted <- subset(gwas_data_lifted, CHRPOS %in% int_chrpos)
    snps_gtex <- snp_attach(subset(snps_gtex,
                                   ind.col = which(snps_gtex$map$CHRPOS %in% int_chrpos)))
    gwas_data_lifted <- gwas_data_lifted[match(snps_gtex$map$CHRPOS,
                                               gwas_data_lifted$CHRPOS), ]
    gwas_data_lifted$Z <- ifelse(gwas_data_lifted$REF == snps_gtex$map$allele2,
                                 gwas_data_lifted$Z,
                                 -1 * gwas_data_lifted$Z)
    
    # Calculate LD matrix
    ld <- as.matrix(snp_cor(snps_gtex$genotypes))
    colnames(ld) <- rownames(ld) <- snps_gtex$map$CHRPOS
    ld_file <- file.path(tempfolder, paste0('LD_', gwas_name, '.txt'))
    write.table(ld, ld_file, sep = '\t', col.names = FALSE, row.names = FALSE)
    
    # Write GWAS Z-scores
    gwas_zfile <- file.path(tempfolder, paste0('gwas_', gwas_name, '.txt'))
    gwas_out <- data.frame(V1 = gwas_data_lifted$CHRPOS, V2 = gwas_data_lifted$Z)
    data.table::fwrite(gwas_out, gwas_zfile, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Process each tissue
    setwd('/rsrch5/scratch/epi/bhattacharya_lab/')
    tissues <- unique(gene_interest$Tissue)
    header_names <- read.table('/rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/nominal_header_0.txt',
                               header = FALSE)
    header_names <- as.character(header_names[1, ])
    header_names = c(header_names,'best_hit')
    header_names[15] = 'std_err'
    
    for (t in tissues) {
      
      cat("    Processing tissue:", t, "\n")
      setwd(file.path('/rsrch5/scratch/epi/bhattacharya_lab/', t))
      dir_out <- file.path('/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare/',
                           t, 'coloc', gwas_name)
      if (any(!file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v45_salmon.txt.gz'))) |
              !file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v38_salmon.txt.gz'))) |
              !file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v27_salmon.txt.gz'))) |
              !file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v45_star.txt.gz'))) |
              !file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v38_star.txt.gz'))) |
              !file.exists(file.path(dir_out, paste0(gene, '_eQTL_coloc_v27_star.txt.gz'))))){
        
        # Read eQTL data for all versions
        v27_star <- data.table::fread(paste0(t, '_v27_star_common_cisQTL_nominal.txt'))
        v27_salmon <- data.table::fread(paste0(t, '_v27_salmon_common_cisQTL_nominal.txt'))
        v38_star <- data.table::fread(paste0(t, '_v38_star_common_cisQTL_nominal.txt'))
        v38_salmon <- data.table::fread(paste0(t, '_v38_salmon_common_cisQTL_nominal.txt'))
        v45_star <- data.table::fread(paste0(t, '_v45_star_common_cisQTL_nominal.txt'))
        v45_salmon <- data.table::fread(paste0(t, '_v45_salmon_common_cisQTL_nominal.txt'))
        colnames(v27_star) <- colnames(v38_star) <- colnames(v45_star) <- header_names
        colnames(v27_salmon) <- colnames(v38_salmon) <- colnames(v45_salmon) <- header_names
        
        v27_star <- subset(v27_star, phe_id == gene)
        v38_star <- subset(v38_star, phe_id == gene)
        v45_star <- subset(v45_star, phe_id == gene)
        v27_salmon <- subset(v27_salmon, phe_id == gene)
        v38_salmon <- subset(v38_salmon, phe_id == gene)
        v45_salmon <- subset(v45_salmon, phe_id == gene)
        
        runColoc = function(qtl_data,name,
                            gwas_name,ldfile,gwas_zfile,
                            gwas_data_lifted){
          
          if (nrow(qtl_data) > 0) {
            qtl_data$CHRPOS <- paste(qtl_data$var_chr, qtl_data$var_from, sep = ':')
            qtl_data <- subset(qtl_data, CHRPOS %in% int_chrpos)
            
            if (nrow(qtl_data) > 0) {
              zfile <- file.path(tempfolder, paste0('eQTL_',name,'_',gwas_name,
                                                    '.txt'))
              qtl_out <- data.frame(V1 = qtl_data$CHRPOS,
                                    V2 = abs(qnorm(qtl_data$nom_pval)) * sign(qtl_data$slope))
              data.table::fwrite(qtl_out, zfile, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
              
              out_ecav <- file.path(tempfolder, paste0('eQTL_',name,'_',
                                                       gwas_name, '_ecav'))
              system(paste(ecaviar,
                           '-o', out_ecav,
                           '-l', ld_file,
                           '-l', ld_file,
                           '-z', gwas_zfile,
                           '-z', zfile,
                           '-c', 2))
              
              if (file.exists(paste0(out_ecav, '_col'))) {
                ecav <- data.table::fread(paste0(out_ecav, '_col'))
                colnames(ecav)[1] <- 'CHRPOS'
                qtl_data <- merge(qtl_data, ecav, by = 'CHRPOS')
                qtl_data <- merge(qtl_data, gwas_data_lifted[, c('CHRPOS', 'P', 'BETA', 'SE', 'Z')], by = 'CHRPOS')
                
                dir_out <- file.path('/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare/',
                                     t, 'coloc', gwas_name)
                dir.create(dir_out, recursive = TRUE)
                
                data.table::fwrite(qtl_data,
                                   file.path(dir_out, paste0(gene, '_eQTL_coloc_',name,'.txt.gz')),
                                   sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
                cat('Written out to',
                    file.path(dir_out, paste0(gene, '_eQTL_coloc_',name,'.txt.gz')))
              }
            }
          }
          
        }
        
        runColoc(qtl_data = v27_star,
                 name = 'v27_star',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        runColoc(qtl_data = v27_salmon,
                 name = 'v27_salmon',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        runColoc(qtl_data = v38_star,
                 name = 'v38_star',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        runColoc(qtl_data = v38_salmon,
                 name = 'v38_salmon',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        runColoc(qtl_data = v45_star,
                 name = 'v45_star',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        runColoc(qtl_data = v45_salmon,
                 name = 'v45_salmon',
                 gwas_name = gwas_name,
                 ldfile = ldfile,
                 gwas_zfile = gwas_zfile,
                 gwas_data_lifted = gwas_data_lifted)
        
      }
    }
    
    # Clean up temp folder for this GWAS
    unlink(tempfolder, recursive = TRUE)
    if (!dir.exists(tempfolder)) {
      cat("    Temp folder", tempfolder, "was deleted!\n")
    }
    
  } # End GWAS loop
  
  cat("Completed gene:", gene, "\n\n")
  
} # End gene loop