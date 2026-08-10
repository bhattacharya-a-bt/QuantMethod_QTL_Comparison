################################################################################
# USER CONFIGURATION - edit these, or set the matching environment variables
################################################################################
# r_lib_path  extra R library path; "" to use the default .libPaths()
# base_path   requantification output root; polished outputs are read from
#             '<base_path>/polished/<annot>/'
# annot       annotation version to count samples for
r_lib_path <- Sys.getenv("R_LIB_PATH",  "")
base_path  <- Sys.getenv("REQUANT_DIR", "/path/to/GTEx_v8/requants")
annot      <- Sys.getenv("ANNOT",       "Ensembl")
################################################################################

if (nchar(r_lib_path) > 0) {
  .libPaths(c(r_lib_path, .libPaths()))
}
library(SummarizedExperiment)

# Polished directory for Ensembl
polished_dir <- file.path(base_path, "polished", annot)

# Identify all tissues with *_featureCounts_gene.RDS
fc_files <- list.files(
  polished_dir,
  pattern = "_featureCounts_gene\\.RDS$",
  full.names = TRUE
)

# Tissue names
tissues <- sort(gsub("_featureCounts_gene\\.RDS$", "", basename(fc_files)))

cat("Found", length(tissues), "tissues in polished directory\n\n")

# Loop and count samples for Ensembl
for (tissue in tissues) {
  cat("------------------------------------------------------------\n")
  cat("TISSUE:", tissue, "\n")
  
  rds_path <- file.path(
    base_path, "polished", annot,
    paste0(tissue, "_featureCounts_gene.RDS")
  )
  
  if (!file.exists(rds_path)) {
    cat("  [Ensembl] missing — skipping\n")
    next
  }
  
  x <- readRDS(rds_path)
  nsamp <- ncol(x[["counts"]])
  
  cat("  samples:", nsamp, "\n\n")
}
