
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
chr <- as.numeric(args[1])

library(data.table)
library(purrr)

annot <- fread("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/anno_selected_genes.txt")
expr_files <- list.files("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/expression")
expr_files <- expr_files[-grep("aggregated",expr_files)]

expr_genes <- expr_files |>
  strsplit("_") |> map_chr(2) |>
  strsplit(".RData") |> map_chr(1) |>
  sub("\\..*$", "", x = _)

dat <- fread("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/results/cis_eqtl/aggregated_per_eqtl_results.txt")
dat <- data.frame(dat)

message("Starting chr ", chr)

bim <- fread(
  paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr,".bim")
) |> data.frame()

dat_chr <- dat[dat$phe_chr == chr,]

dat_chr$true_egene <- dat_chr$gene_id %in% expr_genes
egenes <- unique(dat_chr$gene_id[dat_chr$true_egene])

gene_res <- lapply(seq_along(egenes), function(g){

  gene <- egenes[g]
  if(g %% 50 == 0) message("  gene ", g, "/", length(egenes))

  idx <- which(expr_genes == gene)

  load(paste0(
    "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/expression/",
    expr_files[idx]
  ))

  transcript_info <- do.call(rbind, strsplit(Y$transcript_ids, "\\|"))

  transcript_df <- data.frame(
    transcript_id = sub("\\..*", "", transcript_info[,1]),
    gene_id       = sub("\\..*", "", transcript_info[,2]),
    stringsAsFactors = FALSE
  )
  transcript_df$isoform_h2 <- Y$isoform_h2
  transcript_df$isoform_h2_rank <- order(Y$isoform_h2,decreasing=T)

  B <- Y$B
  colnames(B) <- transcript_df$transcript_id

  snps <- rownames(B)

  snps_pos <- as.integer(map_chr(strsplit(snps,"_"),2))

  bim_sub <- bim[bim$V4 %in% snps_pos,]

  bim_sub$var_id1 <- paste(bim_sub$V1,bim_sub$V4,bim_sub$V5,bim_sub$V6,sep="_")
  bim_sub$var_id2 <- paste(bim_sub$V1,bim_sub$V4,bim_sub$V6,bim_sub$V5,sep="_")

  snps_rs <- sapply(snps, function(x)
    bim_sub$V2[which(bim_sub$var_id1==x | bim_sub$var_id2==x)]
  )

  if(!any(snps_rs == "")){
    rownames(B) <- snps_rs
  } else {
    message("Error matching rsids for ", gene)
  }

  df <- dat_chr[dat_chr$gene_id == gene,]

  res <- df
  res$transcript_id <- vector("list", nrow(res))
  res$beta <- vector("list", nrow(res))
  res$transcript_h2 <- vector("list", nrow(res))
  res$gene_h2 <- Y$gene_h2
  res$gene <- gene
  res$n_transcripts <- ncol(B)
  res$n_qtl_shared <- sum(apply(B,1,FUN=function(x) sum(!x==0)==ncol(B)))
  res$transcript_h2_rank <- vector("list", nrow(res))

  var_results <- lapply(unique(df$var_id), function(var){

    idx_B <- which(rownames(B) == var)
    idx_res <- which(res$var_id == var)

    if(length(idx_B) == 1){

      tmp <- B[idx_B,]
      tmp <- tmp[tmp != 0]

      res$transcript_id[idx_res] <<- list(names(tmp))
      res$beta[idx_res] <<- list(as.numeric(tmp))
      res$transcript_h2[idx_res] <<- list(
        sapply(names(tmp), function(x)
          Y$isoform_h2[transcript_df$transcript_id == x])
      )
      res$transcript_h2_rank[idx_res] <<- list(
        sapply(names(tmp), function(x)
          transcript_df$isoform_h2_rank[transcript_df$transcript_id == x])
      )
    }

  })

  # look up LD between res$var_id and 
  vars <- unique(res$var_id)
  eqtls <- row.names(B)
  all_rsids <- unique(c(vars,eqtls))

  # Define inputs
  plink_prefix <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr)  # without file extensions
  out_prefix <- paste0("/rsrch5/scratch/epi/sthead/tmp/chr",chr,"_",gene)
  snp_file <- tempfile(pattern = "snplist_", fileext = ".txt")
  fwrite(data.table(SNP = all_rsids), snp_file, col.names = FALSE, quote = FALSE)

  # Construct PLINK2 command
  plink_cmd <- paste(
    "plink2",
    "--bfile", plink_prefix,
    "--r2-unphased",
    "--extract", snp_file,
    "--ld-window", 99999,
    "--ld-window-kb", 1000,
    "--ld-window-r2", 0.01,
    "--out", out_prefix
  )

  # Run the command
  system(plink_cmd)

  # Read it in
  vcor_dat <- (fread(paste0(out_prefix,".vcor")))

  # LD between eQTLs and GWAS SNPs
  ld_eqtl_gwas1 <- vcor_dat[
    (ID_A %in% eqtls & ID_B %in% vars)
  ]
  tmp1 <- ld_eqtl_gwas1
  colnames(tmp1) <-c ("true_eqtl_chr","true_eqtl_pos","true_eqtl_rsid","res_chr","res_pos","res_rsid","r2")

  ld_eqtl_gwas2 <- vcor_dat[
    (ID_B %in% eqtls & ID_A %in% vars)
  ]
  tmp2 <- ld_eqtl_gwas2[,c(4,5,6,1,2,3,7)]
  colnames(tmp2) <-c ("true_eqtl_chr","true_eqtl_pos","true_eqtl_rsid","res_chr","res_pos","res_rsid","r2")
  out <- rbind(tmp1,tmp2)

  res$max_ld_with_true_eqtls <- NA
  for(rs in unique(out$res_rsid)){
    max_ld <-  max(out$r2[out$res_rsid==rs])
    res$max_ld_with_true_eqtls[res$var_id==rs] <- max_ld
    rm(max_ld)
  }

  res$shared_eqtl <- (lengths(res$transcript_id) > 0)

  res

})

out <- rbindlist(gene_res)

save(out,file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/results/cis_eqtl/eqtl_res_with_true_beta_chr",chr,".RData"))
