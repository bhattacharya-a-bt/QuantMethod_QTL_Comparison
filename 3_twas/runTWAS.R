#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(bigsnpr)
library(stringr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript script.R <array_index> <gwas_index>")
}

array_index <- as.numeric(args[1])
gwas_index <- as.numeric(args[2])

# Load TWAS model file paths
file_paths <- readRDS('/rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/twas_counts.RDS')
file_paths <- file_paths$File

# Calculate chunk for distributing jobs into 50 arrays
total_files <- length(file_paths)
chunk <- ceiling(total_files / 100)
index_start <- (array_index - 1) * chunk + 1
index_stop <- min(index_start + chunk - 1, total_files)

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

# Select the specific GWAS file for this job
gwas_file <- gwas_files$path[gwas_index]
gwas_name <- gwas_files$name[gwas_index]

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

# Function to extract file information
extract_file_info <- function(file_path) {
  # Split the path by '/'
  path_parts <- str_split(file_path, "/")[[1]]
  
  # Extract tissue (second to last path component)
  tissue <- path_parts[length(path_parts) - 1]
  
  # Extract filename
  filename <- path_parts[length(path_parts)]
  
  # Extract gene (ENSG ID from filename)
  gene <- str_extract(filename, "ENSG\\d+")
  
  # Extract version and method from path
  # Look for pattern like "38_star_common" or "45_salmon_common" in the path
  version_method_part <- str_extract(file_path, "(27|38|45)_(star|salmon)_common")
  
  if (!is.na(version_method_part)) {
    parts <- str_split(version_method_part, "_")[[1]]
    version <- paste0("v", parts[1])
    method <- toupper(parts[2])  # Convert to uppercase
  } else {
    version <- NA
    method <- NA
  }
  
  return(list(
    tissue = tissue,
    gene = gene,
    version = version,
    method = method
  ))
}

# GTEx SNP reference file
gtex_snpfile <- '/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF'

cat("Processing GWAS:", gwas_name, "\n")
cat("Processing array index:", array_index, "covering files", index_start, "to", index_stop, "\n")

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

# Perform liftover on full GWAS data
gwas_data_lifted <- liftover_gwas(gwas_data, chain_file)

# Calculate Z-score if not present
if (!"Z" %in% colnames(gwas_data_lifted)) {
  gwas_data_lifted$Z <- gwas_data_lifted$BETA / gwas_data_lifted$SE
}

# Create CHRPOS column
gwas_data_lifted$CHRPOS <- paste(gwas_data_lifted$CHR, gwas_data_lifted$POS, sep = ':')

# Initialize results data frame
results <- data.frame()

# Loop through each TWAS model file
for (index in index_start:index_stop) {
  
  cat("Processing TWAS model", index, "of", total_files, "\n")
  
  twas_file <- file_paths[index]
  
  # Check if file exists
  if (!file.exists(twas_file)) {
    cat("  TWAS file not found:", twas_file, "\n")
    next
  }
  
  # Read TWAS model
  twas_model <- readRDS(twas_file)
  
  if (nrow(twas_model) == 0) {
    cat("  Empty TWAS model\n")
    next
  }
  
  chr <- twas_model$Chromosome[1]
  min_pos <- min(twas_model$Position)
  max_pos <- max(twas_model$Position)
  
  # Extract file information
  extraction <- extract_file_info(twas_file)
  gene <- extraction$gene
  tissue <- extraction$tissue
  version <- extraction$version
  method <- extraction$method
  
  if (is.na(gene) || is.na(tissue) || is.na(version) || is.na(method)) {
    cat("  Could not extract file information from:", twas_file, "\n")
    next
  }
  
  cat("  Gene:", gene, "Tissue:", tissue, "Version:", version, "Method:", method, "\n")
  
  # Create CHRPOS for TWAS model
  twas_model$CHRPOS <- paste(twas_model$Chromosome, twas_model$Position, sep = ':')
  
  # Create temporary folder
  tempfolder <- paste0('/rsrch5/scratch/epi/abhattacharya3/temp_twas_', array_index, '_', gwas_index, '_', index)
  dir.create(tempfolder, recursive = TRUE, showWarnings = FALSE)
  
  # Use plink2 to extract SNPs for this region
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
  
  # Read genotype data and compute LD
  tryCatch({
    snps <- snp_attach(snp_readBed2(bed_file_path,
                                    backingfile = file.path(tempfolder, paste0(gene, '_', tissue, '.bk'))))
    
    # Create CHRPOS for SNPs
    snps$map$CHRPOS <- paste(snps$map$chromosome, snps$map$physical.pos, sep = ':')
    
    # Compute LD matrix
    ld <- snp_cor(snps$genotypes)
    colnames(ld) <- rownames(ld) <- snps$map$CHRPOS
    
    # Find intersection of GWAS and TWAS model SNPs
    int_chrpos <- intersect(gwas_data_lifted$CHRPOS, twas_model$CHRPOS)
    int_chrpos <- intersect(int_chrpos, snps$map$CHRPOS)
    
    if (length(int_chrpos) == 0) {
      cat("  No overlapping SNPs between GWAS, TWAS model, and genotype data\n")
      unlink(tempfolder, recursive = TRUE)
      next
    }
    
    # Subset data to overlapping SNPs
    gwas_subset <- gwas_data_lifted[gwas_data_lifted$CHRPOS %in% int_chrpos, ]
    twas_subset <- twas_model[twas_model$CHRPOS %in% int_chrpos, ]
    snps_subset_idx <- which(snps$map$CHRPOS %in% int_chrpos)
    ld_subset <- ld[int_chrpos, int_chrpos]
    
    # Order data consistently
    gwas_subset <- gwas_subset[match(int_chrpos, gwas_subset$CHRPOS), ]
    twas_subset <- twas_subset[match(int_chrpos, twas_subset$CHRPOS), ]
    snps_map_subset <- snps$map[snps_subset_idx, ][match(int_chrpos, snps$map[snps_subset_idx, ]$CHRPOS), ]
    
    # Check and flip alleles if necessary
    twas_weights <- twas_subset$Weight
    
    # Compare ALT alleles between GWAS and TWAS model
    # Assuming the reference genotype data matches the TWAS model
    for (i in 1:length(int_chrpos)) {
      if (!is.na(gwas_subset$ALT[i]) && !is.na(snps_map_subset$allele1[i])) {
        if (gwas_subset$ALT[i] != snps_map_subset$allele1[i]) {
          # Flip the weight if alleles don't match
          twas_weights[i] <- -twas_weights[i]
        }
      }
    }
    
    # Compute TWAS Z-statistic
    # Formula: (Weight * GWAS_Z) / sqrt(Weight^T * LD * Weight)
    gwas_z <- gwas_subset$Z
    
    numerator <- sum(twas_weights * gwas_z)
    denominator <- sqrt(as.numeric(t(twas_weights) %*% ld_subset %*% twas_weights))
    
    if (denominator == 0 || is.na(denominator)) {
      cat("  Invalid denominator for TWAS Z-statistic\n")
      twas_z <- NA
    } else {
      twas_z <- numerator / denominator
    }
    
    # Find top GWAS SNP within 1e6 bases
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
    
    # Add result to data frame
    result_row <- data.frame(
      Cancer = gwas_name,
      Gene = gene,
      Tissue = tissue,
      Version = version,
      Method = method,
      TWAS_Z = twas_z,
      Top_GWAS_SNP = top_snp_id,
      Top_GWAS_P = top_snp_pval,
      N_SNPs = length(int_chrpos),
      stringsAsFactors = FALSE
    )
    
    
    # Write out individual result
    dir_out <- file.path('/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare',
                         tissue, 'twas', gwas_name)
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    
    output_file <- file.path(dir_out, paste0(gene, '_TWAS_result_', version, '_', method, '.txt.gz'))
    data.table::fwrite(result_row, output_file, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    cat("  TWAS Z-score:", twas_z, "written to", output_file, "\n")
    
  }, error = function(e) {
    cat("  Error processing:", e$message, "\n")
  })
  
  # Clean up temporary files
  unlink(tempfolder, recursive = TRUE)
}

cat("Completed processing", nrow(results), "TWAS models for", gwas_name, "\n")