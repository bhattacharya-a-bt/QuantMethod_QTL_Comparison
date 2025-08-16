#!/usr/bin/env Rscript
# Load required libraries
require(tximeta)
require(data.table)
require(stringr)
require(edgeR)
require(SummarizedExperiment)

# Get index from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript process_tissue.R <index>")
}
index <- as.integer(args[1])

# Set working directory
setwd('/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess')

# List tissue directories (excluding logs and bam_files_with_sample_attrib.txt)
dirs_all <- list.files()
dirs_all <- dirs_all[!dirs_all %in% c('logs', 'bam_files_with_sample_attrib.txt')]

# Load sample attributes
samp_attrib <- fread('bam_files_with_sample_attrib_rnaseq.txt')

# Get the tissue corresponding to this index
if (index < 1 || index > length(dirs_all)) {
  stop("Invalid index provided.")
}
tissue <- dirs_all[index]
cat("Processing tissue:", tissue, "\n")

# Set working directory to the tissue folder
setwd(file.path('/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess', tissue))

# Subset samples for this tissue
samp_this <- subset(samp_attrib, SMTSD == tissue & !is.na(SMRIN))

# Prepare metadata for tximeta
metadata <- as.data.frame(samp_this)[, c("sample", "SMTSD")]
metadata$names <- unlist(sapply(strsplit(metadata$sample,'-'),function(x) paste0(x[1],'-',x[2])))
metadata$files <- file.path(getwd(), 'salmon_quantifications', 'gencodev45', metadata$sample, 'quant.sf')

# Setup tximeta
setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")
gtf <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"
# makeLinkedTxome(indexDir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45.salmon_index/gencode_v45",
#                 source = "GENCODEv45",
#                 organism = "Homo sapiens",
#                 genome = "GRCh38",
#                 release = "p14",
#                 fasta = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
#                 gtf = gtf,
#                 write = FALSE)

# Read quantifications
se <- tximeta(metadata, type = "salmon", useHub = FALSE, skipSeqinfo = TRUE, txOut = TRUE, dropInfReps = TRUE)

# Summarize to gene-level
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
output_dir <- file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45', tissue)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(se.g, file = file.path(output_dir, paste0(tissue, '_gene.RDS')))
saveRDS(se, file = file.path(output_dir, paste0(tissue, '_transcripts.RDS')))

cat("Finished processing tissue:", tissue, "\n")
cat("Available assays in gene-level object:", names(assays(se.g)), "\n")
cat("Available assays in transcript-level object:", names(assays(se)), "\n")