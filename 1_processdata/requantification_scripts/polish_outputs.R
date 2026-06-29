.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/data/Rlibs/4.3.1", .libPaths()))
library(SummarizedExperiment)

# Polish quantification outputs - find consistent samples and features

ANNOTATION <- c("GENCODE_v27", "GENCODE_v38", "GENCODE_v45", "Ensembl")
base_path <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GTEx_v8/requants"

# Get all RDS files from one annotation directory (they should all have the same tissues)
cat("Identifying all tissues...\n")
files <- list.files(paste0(base_path, "/", ANNOTATION[1]), pattern = "\\.RDS$", full.names = FALSE)

# Extract tissue names by removing the method suffixes
tissues <- unique(gsub("_(featureCounts|kallisto|RSEM|salmon)_(gene|transcripts)\\.RDS$", "", files))

# Sort alphabetically
tissues <- sort(tissues)
cat("Found", length(tissues), "tissues\n\n")

# Loop through each tissue
for (tissue in tissues) {
  cat("================================================================================\n")
  cat("CHECKING TISSUE:", tissue, "\n")
  cat("================================================================================\n\n")
  
  # Check if all input files exist
  all_inputs_exist <- TRUE
  missing_inputs <- c()
  
  for (annot in ANNOTATION) {
    input_files <- c(
      paste0(base_path, "/", annot, "/", tissue, "_featureCounts_gene.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_featureCounts_transcripts.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_kallisto_gene.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_kallisto_transcripts.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_RSEM_gene.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_RSEM_transcripts.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_salmon_gene.RDS"),
      paste0(base_path, "/", annot, "/", tissue, "_salmon_transcripts.RDS")
    )
    
    missing <- input_files[!file.exists(input_files)]
    if (length(missing) > 0) {
      all_inputs_exist <- FALSE
      missing_inputs <- c(missing_inputs, missing)
    }
  }
  
  if (!all_inputs_exist) {
    cat("INPUTS MISSING for", tissue, "- SKIPPING\n")
    cat("Missing files:\n")
    for (mf in missing_inputs) {
      cat("  -", mf, "\n")
    }
    cat("\n")
    next
  }
  
  # Check if all output files already exist
  all_exist <- TRUE
  for (annot in ANNOTATION) {
    files_to_check <- c(
      paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_transcripts.RDS")
    )
    
    if (!all(file.exists(files_to_check))) {
      all_exist <- FALSE
      break
    }
  }
  
  if (all_exist) {
    cat("All output files already exist for", tissue, "- SKIPPING\n\n")
    next
  }
  
  cat("PROCESSING TISSUE:", tissue, "\n\n")
  
  # Check consistent samples and features
  all_keep_samples <- list()
  all_keep_features <- list()
  
  for (annot in ANNOTATION) {
    cat("Working on", annot, "...\n")
    
    cat("  Checking samples and features in featureCounts - gene-level...\n")
    FG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_featureCounts_gene.RDS"))
    FG_samples <- colnames(FG[["counts"]])
    FG_genes <- rownames(FG[["counts"]])
    rm(FG)
    
    cat("  Checking samples and features in featureCounts - transcript-level...\n")
    FT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_featureCounts_transcripts.RDS"))
    FT_samples <- colnames(FT[["counts"]])
    FT_transcripts <- rownames(FT[["counts"]])
    rm(FT)
    
    cat("  Checking samples and features in Kallisto - gene-level...\n")
    KG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_kallisto_gene.RDS"))
    KG_samples <- KG@colData@rownames
    KG_genes <- rownames(KG@assays@data@listData[["counts"]])
    rm(KG)
    
    cat("  Checking samples and features in Kallisto - transcript-level...\n")
    KT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_kallisto_transcripts.RDS"))
    KT_samples <- KT@colData@rownames
    KT_transcripts <- rownames(KT@assays@data@listData[["counts"]])
    rm(KT)
    
    cat("  Checking samples and features in RSEM - gene-level...\n")
    RG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_RSEM_gene.RDS"))
    RG_samples <- colnames(RG[["counts"]])
    RG_genes <- rownames(RG[["counts"]])
    rm(RG)
    
    cat("  Checking samples and features in RSEM - transcript-level...\n")
    RT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_RSEM_transcripts.RDS"))
    RT_samples <- colnames(RT[["counts"]])
    RT_transcripts <- rownames(RT[["counts"]])
    rm(RT)
    
    cat("  Checking samples and features in Salmon - gene-level...\n")
    SG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_salmon_gene.RDS"))
    SG_samples <- SG@colData@rownames
    SG_genes <- rownames(SG@assays@data@listData[["counts"]])
    rm(SG)
    
    cat("  Checking samples and features in Salmon - transcript-level...\n")
    ST <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_salmon_transcripts.RDS"))
    ST_samples <- ST@colData@rownames
    ST_transcripts <- rownames(ST@assays@data@listData[["counts"]])
    rm(ST)
    
    # Find consistent samples across methods within this annotation
    keep_samples <- Reduce(intersect, list(
      FG_samples, FT_samples, KG_samples, KT_samples,
      RG_samples, RT_samples, SG_samples, ST_samples
    ))
    
    # Find consistent features across methods within this annotation
    keep_genes <- Reduce(intersect, list(FG_genes, KG_genes, RG_genes, SG_genes))
    keep_transcripts <- Reduce(intersect, list(FT_transcripts, KT_transcripts, RT_transcripts, ST_transcripts))
    
    all_keep_samples[[annot]] <- keep_samples
    all_keep_features[[annot]] <- list(genes = keep_genes, transcripts = keep_transcripts)
    
    cat("  Samples for", annot, ":", length(keep_samples), "\n")
    cat("  Genes for", annot, ":", length(keep_genes), "\n")
    cat("  Transcripts for", annot, ":", length(keep_transcripts), "\n\n")
  }
  
  # Find intersection of samples across all annotation versions
  keep_samples_final <- Reduce(intersect, all_keep_samples)
  cat("Consistent samples across all annotations:", length(keep_samples_final), "\n\n")
  
  # Check if all output files already exist
  annotations_to_process <- c()
  
  for (annot in ANNOTATION) {
    files_to_check <- c(
      paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_transcripts.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_gene.RDS"),
      paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_transcripts.RDS")
    )
    
    if (!all(file.exists(files_to_check))) {
      annotations_to_process <- c(annotations_to_process, annot)
    }
  }
  
  if (length(annotations_to_process) == 0) {
    cat("All output files already exist for", tissue, "- SKIPPING\n\n")
    next
  }
  
  cat("Annotations to process:", paste(annotations_to_process, collapse=", "), "\n\n")
  
  # Polish
  for (annot in annotations_to_process) {
    cat("Polishing", annot, "...\n")
    
    keep_samples <- keep_samples_final
    keep_genes <- all_keep_features[[annot]][["genes"]]
    keep_transcripts <- all_keep_features[[annot]][["transcripts"]]
    
    cat("  Processing featureCounts - gene-level...\n")
    FG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_featureCounts_gene.RDS"))
    FG[["counts"]] <- FG[["counts"]][keep_genes, keep_samples]
    FG[["TMM_normalized"]] <- FG[["TMM_normalized"]][keep_genes, keep_samples]
    saveRDS(FG, paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_gene.RDS"))
    rm(FG)
    
    cat("  Processing featureCounts - transcript-level...\n")
    FT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_featureCounts_transcripts.RDS"))
    FT[["counts"]] <- FT[["counts"]][keep_transcripts, keep_samples]
    FT[["TMM_normalized"]] <- FT[["TMM_normalized"]][keep_transcripts, keep_samples]
    saveRDS(FT, paste0(base_path, "/polished/", annot, "/", tissue, "_featureCounts_transcripts.RDS"))
    rm(FT)
    
    cat("  Processing Kallisto - gene-level...\n")
    KG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_kallisto_gene.RDS"))
    KG <- KG[keep_genes, keep_samples]
    KG.out <- list(counts=KG@assays@data@listData[["counts"]],
                   TMM_normalized=KG@assays@data@listData[["TMM_normalized"]])
    saveRDS(KG.out, paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_gene.RDS"))
    rm(KG)
    
    cat("  Processing Kallisto - transcript-level...\n")
    KT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_kallisto_transcripts.RDS"))
    KT <- KT[keep_transcripts, keep_samples]
    KT.out <- list(counts=KT@assays@data@listData[["counts"]],
                   TMM_normalized=KT@assays@data@listData[["TMM_normalized"]])
    saveRDS(KT.out, paste0(base_path, "/polished/", annot, "/", tissue, "_kallisto_transcripts.RDS"))
    rm(KT)
    
    cat("  Processing RSEM - gene-level...\n")
    RG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_RSEM_gene.RDS"))
    RG[["counts"]] <- RG[["counts"]][keep_genes, keep_samples]
    RG[["TMM_normalized"]] <- RG[["TMM_normalized"]][keep_genes, keep_samples]
    saveRDS(RG, paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_gene.RDS"))
    rm(RG)
    
    cat("  Processing RSEM - transcript-level...\n")
    RT <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_RSEM_transcripts.RDS"))
    RT[["counts"]] <- RT[["counts"]][keep_transcripts, keep_samples]
    RT[["TMM_normalized"]] <- RT[["TMM_normalized"]][keep_transcripts, keep_samples]
    saveRDS(RT, paste0(base_path, "/polished/", annot, "/", tissue, "_RSEM_transcripts.RDS"))
    rm(RT)
    
    cat("  Processing Salmon - gene-level...\n")
    SG <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_salmon_gene.RDS"))
    SG <- SG[keep_genes, keep_samples]
    SG.out <- list(counts=SG@assays@data@listData[["counts"]],
                   TMM_normalized=SG@assays@data@listData[["TMM_normalized"]])
    saveRDS(SG.out, paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_gene.RDS"))
    rm(SG)
    
    cat("  Processing Salmon - transcript-level...\n")
    ST <- readRDS(paste0(base_path, "/", annot, "/", tissue, "_salmon_transcripts.RDS"))
    ST <- ST[keep_transcripts, keep_samples]
    ST.out <- list(counts=ST@assays@data@listData[["counts"]],
                   TMM_normalized=ST@assays@data@listData[["TMM_normalized"]])
    saveRDS(ST.out, paste0(base_path, "/polished/", annot, "/", tissue, "_salmon_transcripts.RDS"))
    rm(ST)
    
    cat("  Completed", annot, "\n\n")
  }
  
  cat("COMPLETED TISSUE:", tissue, "\n\n")
}

cat("================================================================================\n")
cat("ALL TISSUES PROCESSED SUCCESSFULLY\n")
cat("================================================================================\n")