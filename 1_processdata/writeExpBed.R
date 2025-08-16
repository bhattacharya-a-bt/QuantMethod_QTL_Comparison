# writeExpBed.R
# Process GTEx expression data to create BED files for QTLtools
# Takes an index parameter to select tissue from available directories

print("Starting writeExpBed.R script.")
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

print(paste("Processing index:", index))

# Get list of available tissues
tissue_list <- list.files('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45')
tissue <- tissue_list[index]

print(paste("Selected tissue:", tissue))

biomart_file <- "/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/ensembl_gene_info.RData"

print(paste("Processing tissue:", tissue))

# Install/Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(pacman)
pacman::p_load(dplyr, magrittr, tidyr,
               tibble, ggplot2, janitor,
               DESeq2, SummarizedExperiment,
               bigsnpr, biomaRt, data.table,
               stringr)

# Load human genome GRCh38, coordinates, and biological function (presaved)
load(biomart_file)

# Output directory
out_dir <- paste0("/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/", tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Function to process RDS files (v38 and v45)
process_rds_file <- function(rds_path, version_name) {
  print(paste("Processing", version_name, "file:", rds_path))
  
  if (!file.exists(rds_path)) {
    print(paste("File not found:", rds_path))
    return(NULL)
  }
  
  # Load data
  gene <- readRDS(file = rds_path)
  
  # Get observed genes (remove version number)
  obs_gene <- rownames(gene) %>%
    str_replace(pattern = ".[0-9]+$", replacement = "") %>%
    unique %>%
    data.frame("ensembl_gene_id" = .)
  
  # Join with chromosome and location information
  mapped <- obs_gene %>%
    left_join(gene_info, by = "ensembl_gene_id")
  
  # Filter out high variance/low-mean genes
  # Get expression values (abundance)
  abundance <- assay(gene, "abundance") %>% as.data.frame
  print(paste0('Removed ',sum(duplicated(colnames(abundance))),
               ' samples'))
  abundance = abundance[!duplicated(rownames(abundance)),
                        !duplicated(colnames(abundance))]
  
  # Rule of thumb: discard genes with <= 0.1 TPM for at least 25% of the samples
  stable_genes <- abundance %>%
    filter(rowMeans(. <= 0.1) <= 0.25)  # Keep genes with > 0.1 TPM in > 75% samples
  
  # Add gene expression (abundance) to observed genes
  stable_genes$ensembl_gene_id <- rownames(stable_genes) %>%
    str_replace(pattern = ".[0-9]+$", replacement = "")
  rownames(stable_genes) <- NULL
  stable_genes <- stable_genes %>% relocate(ensembl_gene_id)
  
  # Map expression to gene info
  prebed <- stable_genes %>%
    left_join(mapped, by = "ensembl_gene_id")
  
  # Filter to CHR 1:22
  prebed <- prebed %>%
    filter(chromosome_name %in% as.character(1:22))
  
  # BED format - all genes together (no separation by coding/non-coding)
  bed <- prebed %>%
    dplyr::rename(
      `#Chr` = chromosome_name,
      start = start_position,
      end = end_position,
      pid = ensembl_gene_id
    ) %>%
    mutate(
      gid = pid,
      `#Chr` = as.numeric(`#Chr`)
    ) %>%
    relocate(`#Chr`, start, end, pid, gid, strand) %>%
    arrange(`#Chr`, start, end) %>%
    dplyr::select(-gene_biotype)  # Remove biotype column since we're not separating
  bed[,7:ncol(bed)] = log2(bed[,7:ncol(bed)] + 1)
  
  return(bed)
}

# Function to process v27 bed.gz file
process_v27_file <- function(v27_path) {
  print(paste("Processing v27 file:", v27_path))
  
  if (!file.exists(v27_path)) {
    print(paste("File not found:", v27_path))
    return(NULL)
  }
  
  # Read the compressed bed file
  v27_data <- fread(v27_path)
  
  # The file already has the correct format with #chr, start, end, gene_id
  # We just need to rename columns to match our format and filter to chr 1:22
  bed <- v27_data %>%
    dplyr::rename(
      `#Chr` = `#chr`,
      pid = gene_id
    ) %>%
    mutate(
      gid = pid,
      `#Chr` = as.numeric(str_replace(`#Chr`, "chr", ""))
    ) %>%
    filter(`#Chr` %in% 1:22) %>%
    arrange(`#Chr`, start, end)
  
  first_names = c('#Chr','start','end','pid','gid')
  bed = as.data.frame(bed)
  bed = cbind(bed[,first_names],
              bed[,which(!colnames(bed) %in% first_names)])
  
  return(bed)
}

# Process v38 data
v38_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/", tissue, "/", tissue, "_gene.RDS")
v38_bed <- process_rds_file(v38_path, "v38")

if (!is.null(v38_bed)) {
  v38_output <- paste0(out_dir, "/", tissue, "_v38.bed")
  write.table(v38_bed,
              file = v38_output,
              sep = "\t", 
              quote = FALSE, row.names = FALSE,
              col.names = TRUE)
  
  # Compress and index with bgzip and tabix
  system(paste("bgzip -c", v38_output, ">", paste0(v38_output, ".gz")))
  system(paste("tabix -p bed", paste0(v38_output, ".gz")))
  print(paste("v38 bed file written and indexed:", paste0(v38_output, ".gz")))
}

# Process v45 data
v45_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45/", tissue, "/", tissue, "_gene.RDS")
v45_bed <- process_rds_file(v45_path, "v45")

if (!is.null(v45_bed)) {
  v45_output <- paste0(out_dir, "/", tissue, "_v45.bed")
  write.table(v45_bed,
              file = v45_output,
              sep = "\t", 
              quote = FALSE, row.names = FALSE,
              col.names = TRUE)
  
  # Compress and index with bgzip and tabix
  system(paste("bgzip -c", v45_output, ">", paste0(v45_output, ".gz")))
  system(paste("tabix -p bed", paste0(v45_output, ".gz")))
  print(paste("v45 bed file written and indexed:", paste0(v45_output, ".gz")))
}

# Process v27 data
v27_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gtex_gencodev27/", tissue, ".v8.normalized_expression.bed.gz")
v27_bed <- process_v27_file(v27_path)

if (!is.null(v27_bed)) {
  v27_output <- paste0(out_dir, "/", tissue, "_v27.bed")
  write.table(v27_bed,
              file = v27_output,
              sep = "\t", 
              quote = FALSE, row.names = FALSE,
              col.names = TRUE)
  
  # Compress to .gz
  system(paste("bgzip -c", v27_output, ">", paste0(v27_output, ".gz")))
  system(paste("tabix -p bed", paste0(v27_output, ".gz")))
  print(paste("v27 bed file written and indexed:", paste0(v27_output, ".gz")))
}

print(paste("Script completed successfully for tissue:", tissue))
print(paste("Output files created in:", out_dir))