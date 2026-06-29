#!/usr/bin/env Rscript

########################################################################
# change library to local
########################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

########################################################################
# parse arguments
########################################################################
args <- commandArgs(trailingOnly = TRUE)
pass <- as.numeric(args[1])
param_row_reads <- as.numeric(args[2])
annot <- as.character(args[3])
# countsFromAbundance <- as.character(args[2])

########################################################################
# load dependencies
########################################################################
library(data.table)
library(MASS)
library(pracma)
library(dplyr)
library(SummarizedExperiment)
library(tximeta)
library(tximport)
require(stringr)
require(edgeR)

########################################################################
# load helper functions
########################################################################

### SALMON

files <- list.files(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/salmon/param_row_reads_",param_row_reads,"_",annot),recursive=T,full.names=T)
tmp <- grep("quant.sf",files)
files <- files[tmp]
if(length(files)>498){
  keep <- files
}

coldata_run <- data.frame(
files = keep,
names = basename(dirname(keep))
)

# set cachedir to pass1 regardless of current pass
cache_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass1/files_for_analysis/tximeta")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
}

setTximetaBFC(cache_dir)
# gtf <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"
# makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45",
#                 source = "GENCODEv45",
#                 organism = "Homo sapiens",
#                 genome = "GRCh38",
#                 release = "p14",
#                 fasta = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon_gencode_v45.fasta",
#                 gtf = gtf,
#                 write = T)

# gtf <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.annotation.gtf"
# makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v38.salmon_index",
#                 source = "GENCODEv38",
#                 organism = "Homo sapiens",
#                 genome = "GRCh38",
#                 release = "p14",
#                 fasta = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon_gencode_v38.fasta",
#                 gtf = gtf,
#                 write = T)
# gtf <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/gencode.v27.annotation.gtf"
# makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/salmon/gencode_v27",
#                 source = "GENCODEv27",
#                 organism = "Homo sapiens",
#                 genome = "GRCh38",
#                 release = "p14",
#                 fasta = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v27/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon_gencode_v27.fasta",
#                 gtf = gtf,
#                 write = T)
# gtf <- "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/subset_gtf/gencode.sub.annotation.gtf"
# makeLinkedTxome(indexDir = "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/subset_gtf/index",
#                 source = "GENCODEv38_sub",
#                 organism = "Homo sapiens",
#                 genome = "GRCh38",
#                 release = "p14",
#                 fasta = "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/subset_gtf/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
#                 gtf = gtf,
#                 write = T)


# countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM")
se <- tximeta(coldata_run, skipSeqinfo=T, type="salmon", useHub=F,txOut=T)
se.g <- summarizeToGene(se)

# TMM normalization using edgeR
cat("Performing TMM normalization...\n")

# TMM normalization for gene-level data (se.g)
# Extract count matrix
count_matrix_gene <- assay(se.g, "counts")

# Create DGEList object
dge_gene <- DGEList(counts = count_matrix_gene)

# Calculate normalization factors using TMM
dge_gene <- calcNormFactors(dge_gene, method = "TMM")

# Get TMM-normalized expression values
# Using cpm() which applies the normalization factors
tmm_normalized_gene <- cpm(dge_gene, normalized.lib.sizes = TRUE, log = FALSE)

# Add TMM-normalized expression as a new assay
assay(se.g, "TMM_normalized") <- tmm_normalized_gene

# Optional: also add log2(TMM+1) normalized values
assay(se.g, "TMM_log2") <- log2(tmm_normalized_gene + 1)

cat("TMM normalization completed for gene-level data. Added assays: TMM_normalized, TMM_log2\n")

# TMM normalization for transcript-level data (se)
cat("Performing TMM normalization for transcript-level data...\n")

# Extract count matrix for transcripts
count_matrix_tx <- assay(se, "counts")

# Create DGEList object
dge_tx <- DGEList(counts = count_matrix_tx)

# Calculate normalization factors using TMM
dge_tx <- calcNormFactors(dge_tx, method = "TMM")

# Get TMM-normalized expression values
tmm_normalized_tx <- cpm(dge_tx, normalized.lib.sizes = TRUE, log = FALSE)

# Add TMM-normalized expression as a new assay
assay(se, "TMM_normalized") <- tmm_normalized_tx

# Optional: also add log2(TMM+1) normalized values
assay(se, "TMM_log2") <- log2(tmm_normalized_tx + 1)

cat("TMM normalization completed for transcript-level data. Added assays: TMM_normalized, TMM_log2\n")

# Save outputs
output_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
saveRDS(se.g, file = paste0(output_dir,"/param_row_reads_",param_row_reads,"_salmon_",annot,"_gene.RDS"))
saveRDS(se, file = paste0(output_dir,"/param_row_reads_",param_row_reads,"_salmon_",annot,"_transcripts.RDS"))

cat("Finished processing", "\n")
cat("Available assays in gene-level object:", names(assays(se.g)), "\n")
cat("Available assays in transcript-level object:", names(assays(se)), "\n")



### FEATURECOUNTS

annotation <- annot

output_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/quants")

library(GenomicRanges)
library(edgeR)

se <- readRDS(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/star_alignments/param_row_reads_",param_row_reads,"/",annotation,".RDS"))

counts <- se$counts$counts
annot  <- se$counts$annotation

gr <- GRanges(
  seqnames = annot$Chr,
  ranges = IRanges(start = annot$Start, end = annot$End),
  strand = annot$Strand
)

names(gr) <- annot$GeneID

se.fc <- SummarizedExperiment(
  assays = list(counts = counts),
  rowRanges = gr
)
colnames(se.fc) <- sub("_Aligned.sortedByCoord.out.bam","",colnames(se.fc))

count_matrix_gene <- assay(se.fc, "counts")

dge_gene <- DGEList(counts = count_matrix_gene)
dge_gene <- calcNormFactors(dge_gene, method = "TMM")

tmm_normalized_gene <- cpm(
  dge_gene,
  normalized.lib.sizes = TRUE,
  log = FALSE
)

assay(se.fc, "TMM_normalized") <- tmm_normalized_gene
assay(se.fc, "TMM_log2") <- log2(tmm_normalized_gene + 1)

saveRDS(se.fc, file = paste0(output_dir,"/param_row_reads_",param_row_reads,"_featureCounts_",annotation,"_gene.RDS"))


