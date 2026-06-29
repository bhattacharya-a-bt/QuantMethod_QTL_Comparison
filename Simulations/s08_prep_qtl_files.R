#!/usr/bin/env Rscript

###################################################################
# change library to local
###################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

####################################################################################
# load dependencies
####################################################################################
library(stringr)
library(VariantAnnotation)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(SummarizedExperiment)
library(tximeta)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
pass <- as.numeric(args[1])
param_row_reads <- as.numeric(args[2])
annot <- as.character(args[3])

####################################################################################
# begin code
####################################################################################

# bed file should have the following columns:
# Chromosome ID [string]
# Start genomic position of the phenotype (here the TSS of gene1) [integer, 0-based]
# End genomic position of the phenotype (here the TSS of gene1) [integer, 1-based]
# Phenotype ID (here the exon IDs) [string].
# Phenotype group ID (here the gene IDs, multiple exons belong to the same gene) [string]
# Strand orientation [+/-]

## SALMON

# loop through annotations
for(annot in c("gencode_v38","gencode_v27","gencode_v45","gencode_v38_sub")){
  cat(paste0(annot," \n"))
  se <- readRDS(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants/param_row_reads_",param_row_reads,"_salmon_",annot,"_gene.RDS"))
  rm <- which(colnames(se)=="HG00136")
  if(length(rm)==1){
    se <- se[,-rm]
  }

  tpm_gene <- assay(se, "abundance")

  # Compute proportion of samples > 0.1 TPM
  prop_expr_gene <- rowMeans(tpm_gene > 0.1)
  write.table(prop_expr_gene,quote=F,row.names=T,col.names=F,sep="\t",file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants/param_row_reads_",param_row_reads,"_salmon_",annot,"_prop_samples_ge_01TPM.txt"))

  # filter genes
  se.filt <- se[prop_expr_gene >= 0.25, ]

  gr <- rowRanges(se.filt)
  bed <- data.frame(
    Chr    = as.character(seqnames(gr)),
    start  = start(gr) - 1,              # BED is 0-based
    end    = end(gr),                    # inclusive end
    pid    = names(gr),                  # "pid" – rownames like ENSG...
    gid    = mcols(gr)$gene_id,          # gene IDs
    strand = as.character(strand(gr))
  )

  out_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/bed_files/param_row_reads_",param_row_reads)

  if(!dir.exists(out_dir)){
    dir.create(out_dir,recursive=T)
  }

  # loop through assay names
  for(assay_name in c("abundance","TMM_normalized","TMM_log2")){
    cat(paste0(assay_name," \n"))
    quant <- data.frame(assay(se.filt,assay_name))  
    quant$GeneID <- rownames(quant)

    bed_merged <- merge(x=bed,y=quant,by.x="gid",by.y="GeneID")
    bed_merged <- bed_merged[,c(2,3,4,1,1,6,7:ncol(bed_merged))]
    bed_merged <- bed_merged[bed_merged$Chr %in% paste0("chr",1:22),]
    bed_merged$Chr <- factor(bed_merged$Chr, levels=paste0("chr",1:22),labels=1:22)
    bed_merged <- bed_merged[order(bed_merged$start),]
    bed_merged <- bed_merged[order(bed_merged$Chr),]

    colnames(bed_merged)[1:6] <- c("#Chr","start","end","pid","gid","strand")
    write.table(bed_merged,file=paste0(out_dir,"/salmon_",annot,"_gene_counts_",assay_name,".bed"),sep="\t",quote=F,col.names=T,row.names=F)
  }
}

## FEATURECOUNTS

# loop through annotations
for(annot in c("gencode_v38","gencode_v27","gencode_v45","gencode_v38_sub")){
  cat(paste0(annot," \n"))
  se <- readRDS(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants/param_row_reads_",param_row_reads,"_featureCounts_",annot,"_gene.RDS"))

  rm <- which(colnames(se)=="HG00136")
  if(length(rm)==1){
    se <- se[,-rm]
  }

  counts <- assay(se, "counts")
  gene_length <- width(rowRanges(se))

  rate <- counts / gene_length
  tpm <- t(t(rate) / colSums(rate)) * 1e6

  assay(se, "abundance") <- tpm

  tpm_gene <- assay(se, "abundance")

  # Compute proportion of samples > 0.1 TPM
  prop_expr_gene <- rowMeans(tpm_gene > 0.1)
  write.table(prop_expr_gene,quote=F,row.names=T,col.names=F,sep="\t",file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants/param_row_reads_",param_row_reads,"_featureCounts_",annot,"_prop_samples_ge_01TPM.txt"))

  # filter genes
  se.filt <- se[prop_expr_gene >= 0.25, ]

  gr <- rowRanges(se.filt)
  bed <- data.frame(
    Chr    = as.character(seqnames(gr)),
    start  = start(gr) - 1,              # BED is 0-based
    end    = end(gr),                    # inclusive end
    pid    = names(gr),                  # "pid" – rownames like ENSG...
    gid    = names(gr),          # gene IDs
    strand = as.character(strand(gr))
  )

  out_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/bed_files/param_row_reads_",param_row_reads)

  if(!dir.exists(out_dir)){
    dir.create(out_dir,recursive=T)
  }

  # loop through assay names
  for(assay_name in c("abundance","TMM_normalized","TMM_log2")){
    cat(paste0(assay_name," \n"))
    quant <- data.frame(assay(se.filt,assay_name))  
    quant$GeneID <- rownames(quant)

    bed_merged <- merge(x=bed,y=quant,by.x="gid",by.y="GeneID")
    bed_merged <- bed_merged[,c(2,3,4,1,1,6,7:ncol(bed_merged))]
    bed_merged <- bed_merged[bed_merged$Chr %in% paste0("chr",1:22),]
    bed_merged$Chr <- factor(bed_merged$Chr, levels=paste0("chr",1:22),labels=1:22)
    bed_merged <- bed_merged[order(bed_merged$start),]
    bed_merged <- bed_merged[order(bed_merged$Chr),]

    colnames(bed_merged)[1:6] <- c("#Chr","start","end","pid","gid","strand")
    write.table(bed_merged,file=paste0(out_dir,"/featureCounts_",annot,"_gene_counts_",assay_name,".bed"),sep="\t",quote=F,col.names=T,row.names=F)
  }
}




