################################################################################
# Check GTEx v8 Quantification Output Files
################################################################################
# Purpose: Verify presence of all expected RDS output files across annotations
# Author: stbresnahan
# Date: 2025-11-22
################################################################################

# Base directory
base_dir <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v8/requants"

# Annotation versions to check
annotations <- c("GENCODE_v27", "GENCODE_v38", "GENCODE_v45", "Ensembl")

# Quantification methods
methods <- c("featureCounts", "kallisto", "RSEM", "salmon")

# Levels (gene and transcript)
levels <- c("gene", "transcripts")

################################################################################
# Function to get all unique tissues across annotation directories
################################################################################
get_all_tissues <- function(base_dir, annotations) {
  all_tissues <- c()
  
  for (annot in annotations) {
    annot_dir <- file.path(base_dir, annot)
    if (!dir.exists(annot_dir)) {
      cat("WARNING: Directory does not exist:", annot_dir, "\n")
      next
    }
    
    # Get all RDS files in this annotation directory
    rds_files <- list.files(annot_dir, pattern = "\\.RDS$", full.names = FALSE)
    
    # Extract tissue names from filenames
    # Format: Tissue_Name_method_level.RDS
    for (file in rds_files) {
      # Remove .RDS extension
      base_name <- sub("\\.RDS$", "", file)
      
      # Try to match pattern: tissue_method_level
      # Split by underscore and identify method
      parts <- strsplit(base_name, "_")[[1]]
      
      # Find where the method starts (should be one of our known methods)
      method_idx <- NULL
      for (i in seq_along(parts)) {
        if (parts[i] %in% methods) {
          method_idx <- i
          break
        }
      }
      
      if (!is.null(method_idx) && method_idx > 1) {
        # Tissue name is everything before the method
        tissue <- paste(parts[1:(method_idx-1)], collapse = "_")
        all_tissues <- c(all_tissues, tissue)
      }
    }
  }
  
  # Return unique tissues
  return(unique(all_tissues))
}

################################################################################
# Function to check for missing files
################################################################################
check_missing_files <- function(base_dir, annotations, tissues, methods, levels) {
  
  cat("================================================================================\n")
  cat("GTEx v8 Quantification File Check\n")
  cat("================================================================================\n")
  cat("Base directory:", base_dir, "\n")
  cat("Annotations:", paste(annotations, collapse = ", "), "\n")
  cat("Tissues found:", length(tissues), "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  cat("================================================================================\n\n")
  
  # Track statistics
  total_expected <- 0
  total_found <- 0
  missing_files <- list()
  
  # Check each annotation
  for (annot in annotations) {
    cat("--------------------------------------------------------------------------------\n")
    cat("Annotation:", annot, "\n")
    cat("--------------------------------------------------------------------------------\n")
    
    annot_dir <- file.path(base_dir, annot)
    
    if (!dir.exists(annot_dir)) {
      cat("  ERROR: Directory does not exist!\n\n")
      next
    }
    
    annot_missing <- c()
    annot_expected <- 0
    annot_found <- 0
    
    # Check each tissue
    for (tissue in tissues) {
      tissue_missing <- c()
      
      # Check each method
      for (method in methods) {
        # Check each level
        for (level in levels) {
          # Construct expected filename
          filename <- paste0(tissue, "_", method, "_", level, ".RDS")
          filepath <- file.path(annot_dir, filename)
          
          annot_expected <- annot_expected + 1
          total_expected <- total_expected + 1
          
          # Check if file exists
          if (file.exists(filepath)) {
            annot_found <- annot_found + 1
            total_found <- total_found + 1
          } else {
            tissue_missing <- c(tissue_missing, filename)
            annot_missing <- c(annot_missing, filename)
          }
        }
      }
      
      # Report missing files for this tissue (if any)
      if (length(tissue_missing) > 0) {
        cat("  ", tissue, "- MISSING", length(tissue_missing), "files:\n")
        for (mf in tissue_missing) {
          cat("    -", mf, "\n")
        }
      }
    }
    
    # Summary for this annotation
    cat("\n")
    cat("  Summary for", annot, ":\n")
    cat("    Expected:", annot_expected, "\n")
    cat("    Found:   ", annot_found, "\n")
    cat("    Missing: ", annot_expected - annot_found, "\n")
    
    if (length(annot_missing) == 0) {
      cat("    ✓ All files present!\n")
    }
    
    cat("\n")
    
    # Store missing files for this annotation
    if (length(annot_missing) > 0) {
      missing_files[[annot]] <- annot_missing
    }
  }
  
  # Overall summary
  cat("================================================================================\n")
  cat("OVERALL SUMMARY\n")
  cat("================================================================================\n")
  cat("Total files expected:", total_expected, "\n")
  cat("Total files found:   ", total_found, "\n")
  cat("Total files missing: ", total_expected - total_found, "\n")
  cat("Completion rate:     ", round(100 * total_found / total_expected, 2), "%\n")
  cat("================================================================================\n\n")
  
  # Detailed missing files report
  if (length(missing_files) > 0) {
    cat("DETAILED MISSING FILES REPORT\n")
    cat("================================================================================\n")
    for (annot in names(missing_files)) {
      cat("\n", annot, "(", length(missing_files[[annot]]), "missing ):\n", sep = "")
      for (mf in sort(missing_files[[annot]])) {
        cat("  ", mf, "\n")
      }
    }
    cat("\n")
  } else {
    cat("✓ ALL FILES PRESENT - NO MISSING FILES!\n\n")
  }
  
  return(invisible(missing_files))
}

################################################################################
# Main Execution
################################################################################

# Get all unique tissues
cat("Scanning for tissues...\n")
all_tissues <- get_all_tissues(base_dir, annotations)
cat("Found", length(all_tissues), "unique tissues\n\n")

if (length(all_tissues) > 0) {
  cat("Tissues:\n")
  for (tissue in sort(all_tissues)) {
    cat("  -", tissue, "\n")
  }
  cat("\n")
}

# Check for missing files
missing <- check_missing_files(base_dir, annotations, all_tissues, methods, levels)

# Save results to file
output_file <- file.path(base_dir, "file_check_report.txt")
sink(output_file)

cat("GTEx v8 Quantification File Check Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

if (length(missing) > 0) {
  cat("MISSING FILES BY ANNOTATION:\n\n")
  for (annot in names(missing)) {
    cat(annot, ":\n")
    for (mf in sort(missing[[annot]])) {
      cat("  ", mf, "\n")
    }
    cat("\n")
  }
} else {
  cat("All expected files are present.\n")
}

sink()

cat("Report saved to:", output_file, "\n")

################################################################################
# End of Script
################################################################################