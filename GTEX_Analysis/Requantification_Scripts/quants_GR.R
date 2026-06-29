################################################################################
# GTEx v8 RNA-seq Quantification Pipeline
################################################################################
# Purpose: Process GTEx v8 RNA-seq data across four annotation versions using
#          four different quantification methods
# 
# Author: stbresnahan
# Date: 11-17-2025
# 
# Annotations:
#   - GENCODE v27 (GTEx v8 default)
#   - GENCODE v38
#   - GENCODE v45 (latest GENCODE)
#   - Ensembl v115
#
# Quantification Methods:
#   1. Salmon - Fast, alignment-free transcript quantification
#   2. STAR + featureCounts - Read counting from STAR genome alignments
#   3. STAR + RSEM - Expectation-maximization based quantification from STAR transcriptome alignments
#   4. Kallisto - Pseudoalignment-based transcript quantification
#
# Input requirements:
#   - GTEx v8 sample attributes file
#   - Pre-computed Salmon quantifications (quant.sf files)
#   - STAR-aligned BAM files (for featureCounts)
#   - Pre-computed RSEM results (.genes.results, .isoforms.results)
#   - Pre-computed Kallisto quantifications (abundance.tsv files)
#   - Annotation files (GTF/GFF) for all versions
#   - Salmon/Kallisto indices for each annotation version
#   - SAF annotation files for featureCounts
#
# Output:
#   - Log file at /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/tissue_log.out
#   - RDS files containing gene and transcript-level quantifications
#   - Includes raw counts and TMM-normalized CPM values
#   - Separate outputs for each quantification method and annotation version
#   - featureCounts and RSEM outputs saved as lists with TMM normalization
#   - Salmon and Kallisto outputs saved as SummarizedExperiment objects
################################################################################

# Set library path to personal R library location
.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/data/Rlibs/4.3.1", .libPaths()))

################################################################################
# Load Required Libraries
################################################################################
require(tximeta)                 # For importing Salmon quantifications
require(data.table)              # For efficient data manipulation
require(stringr)                 # For string operations
require(edgeR)                   # For TMM normalization
require(SummarizedExperiment)    # For handling SummarizedExperiment objects
library(Rsubread)                # For featureCounts
library(rtracklayer)             # For importing GTF files
library(pbapply)                 # For parallel processing
library(parallel)                # For parallel processing
library(tximport)                # For importing kallisto quantifications

################################################################################
# Helper Functions for Error Handling
################################################################################

# Safe file existence and size check
check_file_valid <- function(filepath, min_size = 0) {
  if (!file.exists(filepath)) {
    return(list(valid = FALSE, reason = "File does not exist"))
  }
  file_info <- file.info(filepath)
  if (is.na(file_info$size) || file_info$size <= min_size) {
    return(list(valid = FALSE, reason = paste0("File is empty or too small (size: ", file_info$size, " bytes)")))
  }
  return(list(valid = TRUE, reason = "OK"))
}

# Safe RDS save with validation
safe_saveRDS <- function(object, file, object_name = "object") {
  tryCatch({
    if (is.null(object)) {
      return(list(success = FALSE, message = paste0(object_name, " is NULL")))
    }
    if (is.data.frame(object) && nrow(object) == 0) {
      return(list(success = FALSE, message = paste0(object_name, " has 0 rows")))
    }
    if (is.matrix(object) && (nrow(object) == 0 || ncol(object) == 0)) {
      return(list(success = FALSE, message = paste0(object_name, " matrix has 0 dimensions")))
    }
    if (is.list(object) && length(object) == 0) {
      return(list(success = FALSE, message = paste0(object_name, " list is empty")))
    }
    saveRDS(object, file = file)
    # Verify the saved file
    check <- check_file_valid(file, min_size = 100)
    if (!check$valid) {
      return(list(success = FALSE, message = paste0("Saved file validation failed: ", check$reason)))
    }
    return(list(success = TRUE, message = "OK"))
  }, error = function(e) {
    return(list(success = FALSE, message = paste0("Error saving: ", e$message)))
  })
}

# Safe data.table read
safe_fread <- function(filepath, ...) {
  check <- check_file_valid(filepath, min_size = 10)
  if (!check$valid) {
    return(list(data = NULL, success = FALSE, message = check$reason))
  }
  tryCatch({
    dt <- fread(filepath, ...)
    if (nrow(dt) == 0) {
      return(list(data = NULL, success = FALSE, message = "File contains no data rows"))
    }
    return(list(data = dt, success = TRUE, message = "OK"))
  }, error = function(e) {
    return(list(data = NULL, success = FALSE, message = e$message))
  })
}

# Safe TMM normalization with chunked processing
safe_tmm_normalize <- function(count_matrix, name = "matrix", chunk_size = 100) {
  tryCatch({
    if (is.null(count_matrix) || nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
      return(list(normalized = NULL, success = FALSE, message = paste0(name, " is empty or NULL")))
    }
    
    n_samples <- ncol(count_matrix)
    cat("      Total samples:", n_samples, "\n")
    
    # Check for all-zero columns (across entire matrix)
    col_sums <- colSums(count_matrix, na.rm = TRUE)
    if (any(col_sums == 0)) {
      zero_cols <- sum(col_sums == 0)
      warning(paste0("    WARNING: ", zero_cols, " samples have zero total counts in ", name))
    }
    
    # Check for NA/Inf values (across entire matrix)
    if (any(is.na(count_matrix)) || any(is.infinite(count_matrix))) {
      return(list(normalized = NULL, success = FALSE, message = paste0(name, " contains NA or Inf values")))
    }
    
    # If samples <= chunk_size, process normally without chunking
    if (n_samples <= chunk_size) {
      cat("      Processing all", n_samples, "samples at once\n")
      dge <- DGEList(counts = count_matrix)
      dge <- calcNormFactors(dge, method = "TMM")
      tmm_normalized <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
      return(list(normalized = tmm_normalized, success = TRUE, message = "OK"))
    }
    
    # IMPORTANT: Calculate normalization factors on ENTIRE dataset first
    # This ensures consistent normalization across all chunks
    cat("      Calculating TMM normalization factors for all samples...\n")
    dge_full <- DGEList(counts = count_matrix)
    dge_full <- calcNormFactors(dge_full, method = "TMM")
    norm_factors <- dge_full$samples$norm.factors
    lib_sizes <- dge_full$samples$lib.size
    
    # Clear full DGE object to save memory
    rm(dge_full)
    gc()
    
    # Now process in chunks using the pre-calculated normalization factors
    cat("      Processing", n_samples, "samples in chunks of", chunk_size, "\n")
    n_chunks <- ceiling(n_samples / chunk_size)
    
    # Pre-allocate result matrix
    tmm_normalized <- matrix(0, nrow = nrow(count_matrix), ncol = n_samples)
    rownames(tmm_normalized) <- rownames(count_matrix)
    colnames(tmm_normalized) <- colnames(count_matrix)
    
    for (i in 1:n_chunks) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_samples)
      chunk_indices <- start_idx:end_idx
      
      cat("        Chunk", i, "of", n_chunks, "(samples", start_idx, "-", end_idx, ")\n")
      
      # Extract chunk
      count_chunk <- count_matrix[, chunk_indices, drop = FALSE]
      
      # Create DGE object for chunk
      dge_chunk <- DGEList(counts = count_chunk)
      
      # CRITICAL: Apply the pre-calculated normalization factors
      dge_chunk$samples$norm.factors <- norm_factors[chunk_indices]
      dge_chunk$samples$lib.size <- lib_sizes[chunk_indices]
      
      # Calculate CPM using the global normalization factors
      tmm_chunk <- cpm(dge_chunk, normalized.lib.sizes = TRUE, log = FALSE)
      
      # Store results
      tmm_normalized[, chunk_indices] <- tmm_chunk
      
      # Clean up chunk
      rm(count_chunk, dge_chunk, tmm_chunk)
      gc()
    }
    
    cat("      Chunked processing complete\n")
    return(list(normalized = tmm_normalized, success = TRUE, message = "OK"))
    
  }, error = function(e) {
    return(list(normalized = NULL, success = FALSE, message = e$message))
  })
}

################################################################################
# Set Working Parameters
################################################################################
# Get tissue from command line argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No tissue or user specified. Usage: Rscript script.R <tissue> <user name>")
}

tissue <- args[1]
user <- args[2]

cat("================================================================================\n")
cat("GTEx v8 RNA-seq Quantification Pipeline\n")
cat("================================================================================\n")
cat("Target tissue:", tissue, "\n")
cat("User:", user, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

# Working directory for temporary/intermediate files /rsrch5/scratch/epi
wdir <- paste0("/rsrch5/scratch/epi/",user,"/GTEx_v8")
setwd(wdir)
cat("Working directory:", wdir, "\n\n")

################################################################################
# Load and Filter Sample Metadata
################################################################################
cat("Loading GTEx v8 sample metadata...\n")
# Load GTEx v8 sample attributes
metadata_file <- "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx_v8_sample_attributes.txt"
metadata_check <- check_file_valid(metadata_file, min_size = 1000)
if (!metadata_check$valid) {
  stop(paste0("FATAL: Sample metadata file invalid: ", metadata_check$reason))
}

samps <- read.table(
  metadata_file,
  sep = "\t",
  header = TRUE
)

if (nrow(samps) == 0) {
  stop("FATAL: Sample metadata file is empty")
}

# Filter samples to only include specified tissue type
samps.tissue <- samps[samps$tissue == tissue, ]
cat("  Total samples in metadata:", nrow(samps), "\n")
cat("  Samples for", tissue, ":", nrow(samps.tissue), "\n\n")

if (nrow(samps.tissue) == 0) {
  stop(paste0("FATAL: No samples found for tissue '", tissue, "'"))
}

################################################################################
# Main Processing Loop: Process each annotation version
################################################################################
# Process three GENCODE versions (v27, v38, v45) and Ensembl
# For each version, run four quantification methods:
#   1. Salmon (transcript-level, then summarized to gene-level)
#   2. STAR + featureCounts (gene and transcript level)
#   3. STAR + RSEM (gene and transcript level)
#   4. Kallisto (transcript-level, then summarized to gene-level)

# Error and warning logs:
errors <- list()
warnings_log <- list()

cat("================================================================================\n")
cat("Beginning quantification processing\n")
cat("================================================================================\n\n")

for(ANNOTATION in c("GENCODE_v27", "GENCODE_v38", "GENCODE_v45", "Ensembl")) {
  
  cat("--------------------------------------------------------------------------------\n")
  cat("Processing annotation:", ANNOTATION, "\n")
  cat("--------------------------------------------------------------------------------\n")
  
  # Set output directory for this annotation version
  output_dir <- paste0(
    "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/",
    ANNOTATION
  )
  cat("  Output directory:", output_dir, "\n")
  
  # Check output directory exists
  if (!dir.exists(output_dir)) {
    warning_msg <- paste0("Output directory does not exist: ", output_dir)
    cat("  WARNING:", warning_msg, "\n")
    warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "] ", warning_msg)
    tryCatch({
      dir.create(output_dir, recursive = TRUE)
      cat("  Created output directory\n")
    }, error = function(e) {
      errors[[length(errors) + 1]] <- paste0("[", ANNOTATION, "] Failed to create output directory: ", e$message)
      next
    })
  }
  
  # Set annotation-specific paths for Salmon index and GTF annotation
  # Reference genome FASTA (same for all annotations, as all use GRCh38)
  fasta <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta"
  
  # Define annotation-specific parameters
  annotation_params <- list(
    GENCODE_v27 = list(
      index = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/salmon/gencode_v27",
      gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/gencode.v27.annotation.gtf",
      source = "GENCODEv27"
    ),
    GENCODE_v38 = list(
      index = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v38.salmon_index",
      gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.annotation.gtf",
      source = "GENCODEv38"
    ),
    GENCODE_v45 = list(
      index = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45",
      gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf",
      source = "GENCODEv45"
    ),
    Ensembl = list(
      index = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/Ensembl/salmon/Ensembl",
      gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/Ensembl/Ensembl_Homo_sapiens.GRCh38.115.chr.added.gtf",
      source = "Ensembl"
    )
  )
  
  # Extract parameters for current annotation
  index <- annotation_params[[ANNOTATION]]$index
  gtf <- annotation_params[[ANNOTATION]]$gtf
  source <- annotation_params[[ANNOTATION]]$source
  
  # Validate annotation files exist
  gtf_check <- check_file_valid(gtf, min_size = 1000)
  if (!gtf_check$valid) {
    error_msg <- paste0("GTF file invalid for ", ANNOTATION, ": ", gtf_check$reason)
    cat("  ERROR:", error_msg, "\n")
    errors[[length(errors) + 1]] <- error_msg
    cat("  Skipping annotation due to missing GTF\n\n")
    next
  }
  
  cat("  Creating linked transcriptome for tximeta...\n")
  
  # Set tximeta cache location /rsrch5/scratch/epi
  dir.create(paste0("/rsrch5/scratch/epi/", user, "/GTEx_v8/tximeta"), recursive = TRUE)
  setTximetaBFC(paste0("/rsrch5/scratch/epi/",user,"/GTEx_v8/tximeta"))
  
  # Create linked transcriptome for tximeta
  tryCatch({
    makeLinkedTxome(
      indexDir = index,
      source = source,
      organism = "Homo sapiens",
      genome = "GRCh38",
      release = "p14",
      fasta = fasta,
      gtf = gtf,
      write = FALSE
    )
    cat("  Linked transcriptome created successfully\n\n")
  }, error = function(e) {
    warning_msg <- paste0("Failed to create linked transcriptome: ", e$message)
    cat("  WARNING:", warning_msg, "\n")
    warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "] ", warning_msg)
  })
  
  ##############################################################################
  # METHOD 1: Salmon Quantification
  ##############################################################################
  
  # Define output files
  salmon_gene_file <- file.path(output_dir, paste0(tissue, '_salmon_gene.RDS'))
  salmon_tx_file <- file.path(output_dir, paste0(tissue, '_salmon_transcripts.RDS'))
  
  # Only run if both output files don't exist
  if(!file.exists(salmon_gene_file) || !file.exists(salmon_tx_file)) {
    
    cat("  [METHOD 1/4] Processing Salmon quantification...\n")
    
    tryCatch({
      
      # List all Salmon quantification directories
      quant_dir <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/",ANNOTATION,"/", tissue, "/salmon")
      if (!dir.exists(quant_dir)) {
        stop(paste0("Salmon quant directory does not exist: ", quant_dir))
      }
      
      dirs_all <- list.files(quant_dir)
      dirs_all <- dirs_all[dirs_all %in% samps.tissue$SAMPID]
      
      if (length(dirs_all) == 0) {
        stop("No sample directories found matching metadata")
      }
      
      # Prepare metadata data frame for tximeta
      metadata <- samps.tissue[, 1:2]
      metadata$file <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/",
                              ANNOTATION,"/", tissue, "/salmon","/",metadata$SAMPID)
      metadata <- metadata[, c(3, 1, 2)]
      names(metadata) <- c("files", "names", "SMTSD")
      
      # Quality control: Remove samples with missing quant.sf files
      paths <- file.path(metadata$file, "quant.sf")
      missing <- !file.exists(paths)
      if(sum(missing) > 0) {
        metadata <- metadata[!missing, ]
        cat("    Removed", sum(missing), "samples with missing quant.sf files\n")
      }
      
      # Quality control: Check for empty quant.sf files
      if (nrow(metadata) > 0) {
        paths <- file.path(metadata$file, "quant.sf")
        empty_files <- sapply(paths, function(p) {
          fi <- file.info(p)
          is.na(fi$size) || fi$size < 100
        })
        if (sum(empty_files) > 0) {
          metadata <- metadata[!empty_files, ]
          cat("    Removed", sum(empty_files), "samples with empty quant.sf files\n")
        }
      }
      
      # Quality control: Remove samples with missing Salmon metadata
      if (nrow(metadata) > 0) {
        paths <- file.path(metadata$file, "aux_info/meta_info.json")
        missing <- !file.exists(paths)
        if(sum(missing) > 0) {
          metadata <- metadata[!missing, ]
          cat("    Removed", sum(missing), "samples with missing metadata\n")
        }
      }
      
      if(nrow(metadata) == 0) {
        stop("No valid samples remaining after QC")
      }
      
      cat("    Processing", nrow(metadata), "samples\n")
      
      # Update file paths to point directly to quant.sf files
      metadata$files <- paste0(metadata$files, "/quant.sf")
      
      cat("    Importing transcript-level quantifications...\n")
      # Import Salmon quantifications at transcript level
      se <- tximeta(
        metadata,
        type = "salmon",
        useHub = FALSE,
        skipSeqinfo = TRUE,
        txOut = TRUE,
        dropInfReps = TRUE
      )
      
      if (is.null(se) || ncol(se) == 0) {
        stop("tximeta returned empty SummarizedExperiment")
      }
      
      cat("    Summarizing to gene level...\n")
      se.g <- summarizeToGene(se)
      
      if (is.null(se.g) || ncol(se.g) == 0) {
        stop("Gene summarization returned empty result")
      }
      
      cat("    Performing TMM normalization (gene-level)...\n")
      # Gene-level TMM normalization
      count_matrix_gene <- assay(se.g, "counts")
      tmm_result <- safe_tmm_normalize(count_matrix_gene, "gene counts")
      if (!tmm_result$success) {
        warning_msg <- paste0("Gene-level TMM normalization failed: ", tmm_result$message)
        cat("    WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-Salmon] ", warning_msg)
      } else {
        assay(se.g, "TMM_normalized") <- tmm_result$normalized
      }
      
      cat("    Performing TMM normalization (transcript-level)...\n")
      # Transcript-level TMM normalization
      count_matrix_tx <- assay(se, "counts")
      tmm_result_tx <- safe_tmm_normalize(count_matrix_tx, "transcript counts")
      if (!tmm_result_tx$success) {
        warning_msg <- paste0("Transcript-level TMM normalization failed: ", tmm_result_tx$message)
        cat("    WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-Salmon] ", warning_msg)
      } else {
        assay(se, "TMM_normalized") <- tmm_result_tx$normalized
      }
      
      cat("    Saving gene-level results...\n")
      save_result <- safe_saveRDS(se.g, salmon_gene_file, "gene SE")
      if (!save_result$success) {
        stop(paste0("Failed to save gene results: ", save_result$message))
      }
      rm(se.g, count_matrix_gene)
      
      cat("    Saving transcript-level results...\n")
      save_result <- safe_saveRDS(se, salmon_tx_file, "transcript SE")
      if (!save_result$success) {
        stop(paste0("Failed to save transcript results: ", save_result$message))
      }
      rm(se, count_matrix_tx)
      
      cat("    Salmon quantification complete!\n\n")
      rm(tmm_result, tmm_result_tx)
      
    }, error = function(e) {
      error_msg <- paste0("Salmon quantification failed for ", ANNOTATION, ": ", e$message)
      cat("    ERROR:", error_msg, "\n\n")
      errors[[length(errors) + 1]] <<- error_msg
    })
    
  } else {
    cat("  [METHOD 1/4] Skipping Salmon (outputs already exist)\n\n")
  }
  
  
  ##############################################################################
  # METHOD 2: STAR + featureCounts
  ##############################################################################
  
  # Define output files
  fc_gene_file <- file.path(output_dir, paste0(tissue, '_featureCounts_gene.RDS'))
  fc_tx_file <- file.path(output_dir, paste0(tissue, "_featureCounts_transcripts.RDS"))
  
  # Only run if both output files don't exist
  if(!file.exists(fc_gene_file) || !file.exists(fc_tx_file)) {
    
    cat("  [METHOD 2/4] Processing featureCounts...\n")
    
    tryCatch({
      cat("    Reading annotation files...\n")
      
      # Check SAF files exist and are valid
      gene_saf <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/fc_annots/",ANNOTATION,"_genes.saf")
      tx_saf <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/fc_annots/",ANNOTATION,"_transcripts.saf")
      
      gene_saf_check <- check_file_valid(gene_saf, min_size = 100)
      if (!gene_saf_check$valid) {
        stop(paste0("Gene SAF file invalid: ", gene_saf_check$reason))
      }
      
      tx_saf_check <- check_file_valid(tx_saf, min_size = 100)
      if (!tx_saf_check$valid) {
        stop(paste0("Transcript SAF file invalid: ", tx_saf_check$reason))
      }
      
      # Read gene annotations
      genes <- read.table(gene_saf, sep="\t", header=T)
      transcripts <- read.table(tx_saf, sep="\t", header=T)
      
      if (nrow(genes) == 0) {
        stop("Gene SAF file is empty")
      }
      if (nrow(transcripts) == 0) {
        stop("Transcript SAF file is empty")
      }
      
      # Prepare count matrices
      fC.g <- data.frame(GeneID=unique(genes$GeneID))
      fC.tx <- data.frame(GeneID=unique(transcripts$GeneID))
      
      # Path to the featureCounts directory
      fc_dir <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/",
                       ANNOTATION,"/",tissue,"/featureCounts")
      
      if (!dir.exists(fc_dir)) {
        stop(paste0("featureCounts directory does not exist: ", fc_dir))
      }
      
      cat("    Processing gene-level counts...\n")
      # All gene count files
      files <- list.files(fc_dir, pattern = "_gene.txt$", full.names = TRUE)
      
      if(length(files) == 0) {
        stop("No gene count files found")
      }
      
      # Validate files are not empty
      valid_files <- sapply(files, function(f) {
        fi <- file.info(f)
        !is.na(fi$size) && fi$size > 100
      })
      if (sum(!valid_files) > 0) {
        cat("      Removed", sum(!valid_files), "empty gene count files\n")
        files <- files[valid_files]
      }
      
      if (length(files) == 0) {
        stop("No valid gene count files after filtering")
      }
      
      cat("      Found", length(files), "valid gene count files\n")
      
      # Load base table
      fC.g <- as.data.table(fC.g)
      fC.g[, GeneID := as.character(GeneID)]
      
      # Process gene counts files in parallel
      n_cores <- min(8, parallel::detectCores() - 1)
      cat("      Processing files in parallel (", n_cores, "cores)...\n")
      cl <- makeCluster(n_cores)
      clusterEvalQ(cl, .libPaths(c("/rsrch5/home/epi/bhattacharya_lab/data/Rlibs/4.3.1", .libPaths())))
      clusterEvalQ(cl, library(data.table))
      
      # Process files in parallel with error handling
      res_list <- pblapply(files, function(f) {
        tryCatch({
          fname <- basename(f)
          sample_id <- sub("_gene\\.txt$", "", fname)
          dt <- fread(f, sep = "\t", skip = 1)
          if (nrow(dt) == 0) return(NULL)
          count_col <- names(dt)[ncol(dt)]
          out <- dt[, .(GeneID = Geneid)]
          out[[sample_id]] <- dt[[count_col]]
          out
        }, error = function(e) {
          return(NULL)
        })
      }, cl = cl)
      
      stopCluster(cl)
      
      # Remove NULL results
      res_list <- res_list[!sapply(res_list, is.null)]
      
      if (length(res_list) == 0) {
        stop("All gene count files failed to process")
      }
      
      cat("      Successfully processed", length(res_list), "files\n")
      cat("      Merging count matrices...\n")
      
      # Merge all into one data.table
      merged_counts <- res_list[[1]]
      setkey(merged_counts, GeneID)
      for (i in 2:length(res_list)) {
        setkey(res_list[[i]], GeneID)
        merged_counts <- res_list[[i]][merged_counts]
      }
      
      fC.g <- merge(fC.g, merged_counts, by = "GeneID", all.x = TRUE)
      setnames(fC.g, 1, "gene_id")
      
      cat("      Performing TMM normalization...\n")
      gene_ids <- fC.g$gene_id
      count_matrix_gene <- as.matrix(fC.g[, -1])
      rownames(count_matrix_gene) <- gene_ids
      
      # Replace NA with 0 for TMM
      count_matrix_gene[is.na(count_matrix_gene)] <- 0
      
      tmm_result <- safe_tmm_normalize(count_matrix_gene, "featureCounts gene")
      if (!tmm_result$success) {
        warning_msg <- paste0("Gene TMM normalization failed: ", tmm_result$message)
        cat("      WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-featureCounts] ", warning_msg)
        tmm_normalized_gene <- NULL
      } else {
        tmm_normalized_gene <- tmm_result$normalized
      }
      
      fC.g_list <- list(
        counts = count_matrix_gene,
        TMM_normalized = tmm_normalized_gene
      )
      
      cat("      Saving gene-level results...\n")
      save_result <- safe_saveRDS(fC.g_list, fc_gene_file, "featureCounts gene list")
      if (!save_result$success) {
        stop(paste0("Failed to save gene results: ", save_result$message))
      }
      rm(fC.g, fC.g_list, genes, merged_counts, res_list, count_matrix_gene, tmm_normalized_gene)
      
      
      cat("    Processing transcript-level counts...\n")
      files <- list.files(fc_dir, pattern = "_transcript.txt$", full.names = TRUE)
      
      if (length(files) == 0) {
        stop("No transcript count files found")
      }
      
      # Validate files
      valid_files <- sapply(files, function(f) {
        fi <- file.info(f)
        !is.na(fi$size) && fi$size > 100
      })
      if (sum(!valid_files) > 0) {
        cat("      Removed", sum(!valid_files), "empty transcript count files\n")
        files <- files[valid_files]
      }
      
      if (length(files) == 0) {
        stop("No valid transcript count files after filtering")
      }
      
      cat("      Found", length(files), "valid transcript count files\n")
      
      fC.tx <- as.data.table(fC.tx)
      fC.tx[, GeneID := as.character(GeneID)]
      
      n_cores <- min(8, parallel::detectCores() - 1)
      cat("      Processing files in parallel (", n_cores, "cores)...\n")
      cl <- makeCluster(n_cores)
      clusterEvalQ(cl, .libPaths(c("/rsrch5/home/epi/bhattacharya_lab/data/Rlibs/4.3.1", .libPaths())))
      clusterEvalQ(cl, library(data.table))
      
      res_list <- pblapply(files, function(f) {
        tryCatch({
          fname <- basename(f)
          sample_id <- sub("_transcript\\.txt$", "", fname)
          dt <- fread(f, sep = "\t", skip = 1)
          if (nrow(dt) == 0) return(NULL)
          count_col <- names(dt)[ncol(dt)]
          out <- dt[, .(GeneID = Geneid)]
          out[[sample_id]] <- dt[[count_col]]
          out
        }, error = function(e) {
          return(NULL)
        })
      }, cl = cl)
      
      stopCluster(cl)
      
      res_list <- res_list[!sapply(res_list, is.null)]
      
      if (length(res_list) == 0) {
        stop("All transcript count files failed to process")
      }
      
      cat("      Successfully processed", length(res_list), "files\n")
      cat("      Merging count matrices...\n")
      
      merged_counts <- res_list[[1]]
      setkey(merged_counts, GeneID)
      for (i in 2:length(res_list)) {
        setkey(res_list[[i]], GeneID)
        merged_counts <- res_list[[i]][merged_counts]
      }
      
      fC.tx <- merge(fC.tx, merged_counts, by = "GeneID", all.x = TRUE)
      setnames(fC.tx, 1, "transcript_id")
      
      cat("      Performing TMM normalization...\n")
      tx_ids <- fC.tx$transcript_id
      count_matrix_tx <- as.matrix(fC.tx[, -1])
      rownames(count_matrix_tx) <- tx_ids
      count_matrix_tx[is.na(count_matrix_tx)] <- 0
      
      tmm_result <- safe_tmm_normalize(count_matrix_tx, "featureCounts transcript")
      if (!tmm_result$success) {
        warning_msg <- paste0("Transcript TMM normalization failed: ", tmm_result$message)
        cat("      WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-featureCounts] ", warning_msg)
        tmm_normalized_tx <- NULL
      } else {
        tmm_normalized_tx <- tmm_result$normalized
      }
      
      fC.tx_list <- list(
        counts = count_matrix_tx,
        TMM_normalized = tmm_normalized_tx
      )
      
      cat("      Saving transcript-level results...\n")
      save_result <- safe_saveRDS(fC.tx_list, fc_tx_file, "featureCounts transcript list")
      if (!save_result$success) {
        stop(paste0("Failed to save transcript results: ", save_result$message))
      }
      rm(fC.tx, fC.tx_list, transcripts, merged_counts, res_list, count_matrix_tx, tmm_normalized_tx)
      
      cat("    featureCounts complete!\n\n")
      
    }, error = function(e) {
      error_msg <- paste0("featureCounts failed for ", ANNOTATION, ": ", e$message)
      cat("    ERROR:", error_msg, "\n\n")
      errors[[length(errors) + 1]] <<- error_msg
    })
    
  } else {
    cat("  [METHOD 2/4] Skipping featureCounts (outputs already exist)\n\n")
  }
  
  ##############################################################################
  # METHOD 3: STAR + RSEM
  ##############################################################################
  
  rsem_input_dir <- file.path(
    "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants",
    ANNOTATION, tissue,
    "RSEM"
  )
  
  rsem_gene_file <- file.path(output_dir, paste0(tissue, "_RSEM_gene.RDS"))
  rsem_tx_file <- file.path(output_dir, paste0(tissue, "_RSEM_transcripts.RDS"))
  
  if(!file.exists(rsem_gene_file) || !file.exists(rsem_tx_file)) {
    
    cat("  [METHOD 3/4] Processing RSEM...\n")
    
    tryCatch({
      if (!dir.exists(rsem_input_dir)) {
        stop(paste0("RSEM directory does not exist: ", rsem_input_dir))
      }
      
      gene_files <- list.files(rsem_input_dir, pattern = "\\.genes.results$", full.names = TRUE)
      tx_files <- list.files(rsem_input_dir, pattern = "\\.isoforms.results$", full.names = TRUE)
      
      if (length(gene_files) == 0) {
        stop("No RSEM gene result files found")
      }
      
      # Validate files are not empty
      valid_gene <- sapply(gene_files, function(f) {
        fi <- file.info(f)
        !is.na(fi$size) && fi$size > 100
      })
      if (sum(!valid_gene) > 0) {
        cat("    Removed", sum(!valid_gene), "empty gene result files\n")
        gene_files <- gene_files[valid_gene]
      }
      
      valid_tx <- sapply(tx_files, function(f) {
        fi <- file.info(f)
        !is.na(fi$size) && fi$size > 100
      })
      if (sum(!valid_tx) > 0) {
        cat("    Removed", sum(!valid_tx), "empty isoform result files\n")
        tx_files <- tx_files[valid_tx]
      }
      
      if (length(gene_files) == 0) {
        stop("No valid RSEM gene files after filtering")
      }
      
      cat("    Found", length(gene_files), "valid gene result files\n")
      cat("    Found", length(tx_files), "valid isoform result files\n")
      
      # Validate column names in gene files
      cat("    Validating gene file column names...\n")
      required_gene_cols <- c("gene_id", "expected_count", "TPM", "effective_length")
      valid_gene_cols <- sapply(gene_files, function(f) {
        tryCatch({
          header <- read.table(f, header = TRUE, nrows = 1)
          all(required_gene_cols %in% names(header))
        }, error = function(e) FALSE)
      })
      
      if (sum(!valid_gene_cols) > 0) {
        cat("    Removed", sum(!valid_gene_cols), "gene files with invalid columns\n")
        gene_files <- gene_files[valid_gene_cols]
      }
      
      if (length(gene_files) == 0) {
        stop("No valid RSEM gene files after column validation")
      }
      
      # Validate column names in transcript files
      if (length(tx_files) > 0) {
        cat("    Validating transcript file column names...\n")
        required_tx_cols <- c("transcript_id", "expected_count", "TPM", "effective_length")
        valid_tx_cols <- sapply(tx_files, function(f) {
          tryCatch({
            header <- read.table(f, header = TRUE, nrows = 1)
            all(required_tx_cols %in% names(header))
          }, error = function(e) FALSE)
        })
        
        if (sum(!valid_tx_cols) > 0) {
          cat("    Removed", sum(!valid_tx_cols), "transcript files with invalid columns\n")
          tx_files <- tx_files[valid_tx_cols]
        }
      }
      
      # Validate consistent gene sets across samples
      cat("    Validating consistent gene sets across samples...\n")
      
      num_cores <- round(detectCores()/2)
      
      gene_ids_list <- mclapply(gene_files, function(f) {
        tryCatch({
          df <- read.table(f, header = TRUE, nrows = -1)
          df$gene_id
        }, error = function(e) character(0))
      }, mc.cores = num_cores)
      
      # Check all files have same number of genes
      gene_counts <- sapply(gene_ids_list, length)
      if (length(unique(gene_counts)) > 1) {
        cat("    WARNING: Files have different numbers of genes:\n")
        print(table(gene_counts))
        
        # Keep only files with the most common gene count
        most_common_count <- as.numeric(names(sort(table(gene_counts), decreasing = TRUE)[1]))
        valid_files <- gene_counts == most_common_count
        
        cat("    Removed", sum(!valid_files), "files with inconsistent gene counts\n")
        gene_files <- gene_files[valid_files]
        
        if (length(gene_files) == 0) {
          stop("No valid files remaining after gene count validation")
        }
      }
      
      # Validate consistent transcript sets if applicable
      if (length(tx_files) > 0) {
        cat("    Validating consistent transcript sets across samples...\n")
        num_cores <- round(detectCores()/2)
        tx_ids_list <- mclapply(tx_files, function(f) {
          tryCatch({
            df <- read.table(f, header = TRUE, nrows = -1)
            df$transcript_id
          }, error = function(e) character(0))
        }, mc.cores = num_cores)
        
        # Check all files have same number of transcripts
        tx_counts <- sapply(tx_ids_list, length)
        if (length(unique(tx_counts)) > 1) {
          cat("    WARNING: Files have different numbers of transcripts:\n")
          print(table(tx_counts))
          
          # Keep only files with the most common transcript count
          most_common_count <- as.numeric(names(sort(table(tx_counts), decreasing = TRUE)[1]))
          valid_files <- tx_counts == most_common_count
          
          cat("    Removed", sum(!valid_files), "files with inconsistent transcript counts\n")
          tx_files <- tx_files[valid_files]
        }
      }
      
      # Extract sample names AFTER all validation
      gene_sample_names <- gsub("\\.genes.results$", "", basename(gene_files))
      names(gene_files) <- gene_sample_names
      
      if (length(tx_files) > 0) {
        tx_sample_names <- gsub("\\.isoforms.results$", "", basename(tx_files))
        names(tx_files) <- tx_sample_names
      }
      
      cat("    Importing gene-level quantifications...\n")
      gene_tximport <- tximport(files = gene_files,
                                type = "rsem",
                                txIn = FALSE,
                                txOut = FALSE)
      
      if (is.null(gene_tximport$counts) || ncol(gene_tximport$counts) == 0) {
        stop("Gene-level import returned empty counts")
      }
      
      cat("    Performing TMM normalization (gene-level)...\n")
      count_matrix_gene <- gene_tximport$counts
      
      tmm_result <- safe_tmm_normalize(count_matrix_gene, "RSEM gene")
      if (!tmm_result$success) {
        warning_msg <- paste0("Gene TMM normalization failed: ", tmm_result$message)
        cat("    WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-RSEM] ", warning_msg)
        tmm_normalized_gene <- NULL
      } else {
        tmm_normalized_gene <- tmm_result$normalized
      }
      
      gene_output <- list(
        counts = count_matrix_gene,
        TMM_normalized = tmm_normalized_gene
      )
      
      cat("    Saving gene-level results...\n")
      save_result <- safe_saveRDS(gene_output, rsem_gene_file, "RSEM gene list")
      if (!save_result$success) {
        stop(paste0("Failed to save gene results: ", save_result$message))
      }
      rm(gene_tximport, count_matrix_gene, tmm_normalized_gene, gene_output)
      
      if (length(tx_files) > 0) {
        cat("    Importing transcript-level quantifications...\n")
        tx_tximport <- tximport(files = tx_files,
                                type = "rsem",
                                txIn = TRUE,
                                txOut = TRUE)
        
        if (is.null(tx_tximport$counts) || ncol(tx_tximport$counts) == 0) {
          stop("Transcript-level import returned empty counts")
        }
        
        cat("    Performing TMM normalization (transcript-level)...\n")
        count_matrix_tx <- tx_tximport$counts
        
        tmm_result <- safe_tmm_normalize(count_matrix_tx, "RSEM transcript")
        if (!tmm_result$success) {
          warning_msg <- paste0("Transcript TMM normalization failed: ", tmm_result$message)
          cat("    WARNING:", warning_msg, "\n")
          warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-RSEM] ", warning_msg)
          tmm_normalized_tx <- NULL
        } else {
          tmm_normalized_tx <- tmm_result$normalized
        }
        
        tx_output <- list(
          counts = count_matrix_tx,
          TMM_normalized = tmm_normalized_tx
        )
        
        cat("    Saving transcript-level results...\n")
        save_result <- safe_saveRDS(tx_output, rsem_tx_file, "RSEM transcript list")
        if (!save_result$success) {
          stop(paste0("Failed to save transcript results: ", save_result$message))
        }
        rm(tx_tximport, count_matrix_tx, tmm_normalized_tx, tx_output)
      }
      
      cat("    RSEM complete!\n\n")
      rm(gene_files, tx_files, gene_sample_names, tx_sample_names, 
         gene_ids_list, tx_ids_list, gene_tximport, tx_tximport,
         count_matrix_gene, count_matrix_tx, tmm_result,
         gene_output, tx_output, tmm_normalized_gene, tmm_normalized_tx)
      
    }, error = function(e) {
      error_msg <- paste0("RSEM failed for ", ANNOTATION, ": ", e$message)
      cat("    ERROR:", error_msg, "\n\n")
      errors[[length(errors) + 1]] <<- error_msg
    })
    
  } else {
    cat("  [METHOD 3/4] Skipping RSEM (outputs already exist)\n\n")
  }
  
  
  ##############################################################################
  # METHOD 4: Kallisto Quantification
  ##############################################################################
  
  kallisto_gene_file <- file.path(output_dir, paste0(tissue, "_kallisto_gene.RDS"))
  kallisto_tx_file <- file.path(output_dir, paste0(tissue, "_kallisto_transcripts.RDS"))
  
  if(!file.exists(kallisto_gene_file) || !file.exists(kallisto_tx_file)) {
    
    cat("  [METHOD 4/4] Processing Kallisto quantification...\n")
    
    tryCatch({
      cat("    Reading transcript-to-gene mapping...\n")
      tx2g_file_map <- c(
        GENCODE_v27="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/GENCODE_v27_tx2gene.csv",
        GENCODE_v38="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/GENCODE_v38_tx2gene.csv",
        GENCODE_v45="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/GENCODE_v45_tx2gene.csv",
        Ensembl="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/Ensembl/Ensembl_tx2gene.csv"
      )
      tx2gene_file <- tx2g_file_map[ANNOTATION]
      
      tx2gene_check <- check_file_valid(tx2gene_file, min_size = 100)
      if (!tx2gene_check$valid) {
        stop(paste0("tx2gene file invalid: ", tx2gene_check$reason))
      }
      
      tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)
      if (nrow(tx2gene) == 0) {
        stop("tx2gene file is empty")
      }
      
      # /rsrch5/scratch/epi
      kallisto_dir <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/",ANNOTATION,"/", tissue, "/kallisto")
      if (!dir.exists(kallisto_dir)) {
        stop(paste0("Kallisto directory does not exist: ", kallisto_dir))
      }
      
      dirs_all <- list.files(kallisto_dir)
      dirs_all <- dirs_all[dirs_all %in% samps.tissue$SAMPID]
      
      if (length(dirs_all) == 0) {
        stop("No sample directories found matching metadata")
      }
      
      files <- file.path("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants",
                         ANNOTATION, tissue, "kallisto", dirs_all, "abundance.tsv")
      names(files) <- dirs_all
      
      # Remove missing abundance.tsv
      missing_abund <- !file.exists(files)
      if (sum(missing_abund) > 0) {
        cat("    Removed", sum(missing_abund), "samples with missing abundance.tsv files\n")
        files <- files[!missing_abund]
      }
      
      # Check for empty files
      if (length(files) > 0) {
        empty_files <- sapply(files, function(f) {
          fi <- file.info(f)
          is.na(fi$size) || fi$size < 100
        })
        if (sum(empty_files) > 0) {
          cat("    Removed", sum(empty_files), "samples with empty abundance.tsv files\n")
          files <- files[!empty_files]
        }
      }
      
      # Check run_info.json
      if (length(files) > 0) {
        runinfo <- file.path(dirname(files), "run_info.json")
        names(runinfo) <- names(files)
        missing_runinfo <- !file.exists(runinfo)
        if (sum(missing_runinfo) > 0) {
          cat("    Removed", sum(missing_runinfo), "samples with missing run_info.json files\n")
          files <- files[!missing_runinfo]
        }
      }
      
      if (length(files) == 0) {
        stop("No valid samples remaining after QC")
      }
      
      cat("    Processing", length(files), "samples\n")
      
      cat("    Importing transcript-level quantifications...\n")
      txi <- tximport(files,
                      type = "kallisto",
                      txOut = TRUE,
                      tx2gene = tx2gene,
                      ignoreTxVersion = TRUE)
      
      if (is.null(txi$counts) || ncol(txi$counts) == 0) {
        stop("Kallisto import returned empty counts")
      }
      
      cat("    Creating SummarizedExperiment object...\n")
      se <- SummarizedExperiment(
        assays = list(counts = txi$counts,
                      abundance = txi$abundance,
                      length = txi$length),
        colData = samps.tissue[match(names(files), samps.tissue$SAMPID), ]
      )
      
      cat("    Summarizing to gene level...\n")
      txi.g <- tximport::summarizeToGene(txi, tx2gene = tx2gene)
      se.g <- SummarizedExperiment(
        assays = list(
          counts = txi.g$counts,
          abundance = txi.g$abundance,
          length = txi.g$length
        ),
        colData = samps.tissue[match(names(files), samps.tissue$SAMPID), ]
      )
      
      cat("    Performing TMM normalization (gene-level)...\n")
      count_matrix_gene <- assay(se.g, "counts")
      tmm_result <- safe_tmm_normalize(count_matrix_gene, "Kallisto gene")
      if (!tmm_result$success) {
        warning_msg <- paste0("Gene TMM normalization failed: ", tmm_result$message)
        cat("    WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-Kallisto] ", warning_msg)
      } else {
        assay(se.g, "TMM_normalized") <- tmm_result$normalized
      }
      
      cat("    Performing TMM normalization (transcript-level)...\n")
      count_matrix_tx <- assay(se, "counts")
      tmm_result <- safe_tmm_normalize(count_matrix_tx, "Kallisto transcript")
      if (!tmm_result$success) {
        warning_msg <- paste0("Transcript TMM normalization failed: ", tmm_result$message)
        cat("    WARNING:", warning_msg, "\n")
        warnings_log[[length(warnings_log) + 1]] <- paste0("[", ANNOTATION, "-Kallisto] ", warning_msg)
      } else {
        assay(se, "TMM_normalized") <- tmm_result$normalized
      }
      
      cat("    Saving gene-level results...\n")
      save_result <- safe_saveRDS(se.g, kallisto_gene_file, "Kallisto gene SE")
      if (!save_result$success) {
        stop(paste0("Failed to save gene results: ", save_result$message))
      }
      rm(se.g, count_matrix_gene)
      
      cat("    Saving transcript-level results...\n")
      save_result <- safe_saveRDS(se, kallisto_tx_file, "Kallisto transcript SE")
      if (!save_result$success) {
        stop(paste0("Failed to save transcript results: ", save_result$message))
      }
      rm(se, count_matrix_tx)
      
      cat("    Kallisto quantification complete!\n\n")
      rm(tx2gene, files, dirs_all, txi, txi.g, se, se.g, 
         count_matrix_gene, count_matrix_tx, tmm_result,
         missing_abund, empty_files, runinfo, missing_runinfo)
      
    }, error = function(e) {
      error_msg <- paste0("Kallisto failed for ", ANNOTATION, ": ", e$message)
      cat("    ERROR:", error_msg, "\n\n")
      errors[[length(errors) + 1]] <<- error_msg
    })
    
  } else {
    cat("  [METHOD 4/4] Skipping Kallisto (outputs already exist)\n\n")
  }
  
  cat("Completed processing for", ANNOTATION, "\n")
  cat("--------------------------------------------------------------------------------\n\n")
  
}

################################################################################
# Write Log File
################################################################################

log_messages <- c()

if (length(errors) > 0) {
  log_messages <- c("Pipeline completed with ERRORS:", "", "ERRORS:", unlist(errors))
} else {
  log_messages <- c("Pipeline completed successfully!")
}

if (length(warnings_log) > 0) {
  log_messages <- c(log_messages, "", "WARNINGS:", unlist(warnings_log))
}

log_messages <- c(log_messages, "", paste0("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

df <- data.frame(message = log_messages, stringsAsFactors = FALSE)

log_file <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/", tissue, "_combine.log")
write.table(
  df,
  file = log_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("================================================================================\n")
cat("Pipeline completed\n")
if (length(errors) > 0) {
  cat("  ERRORS:", length(errors), "\n")
}
if (length(warnings_log) > 0) {
  cat("  WARNINGS:", length(warnings_log), "\n")
}
cat("See", log_file, "for details\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n")

################################################################################
# End of Script
################################################################################