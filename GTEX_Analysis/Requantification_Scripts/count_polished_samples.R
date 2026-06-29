.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/data/Rlibs/4.3.1", .libPaths()))
library(SummarizedExperiment)

base_path <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v8/requants"
annot     <- "Ensembl"

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
