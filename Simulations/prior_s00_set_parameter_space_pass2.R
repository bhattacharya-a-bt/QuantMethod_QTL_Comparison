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
library(VariantAnnotation)
library(dplyr)
library(rtracklayer)
library(Rsubread)

####################################################################################
# set static parameters
####################################################################################
pass <- 2
min_num_tx <- 1
num_genes_selected <- "all" # can be "all" or number < ~12k
seed <- 123
n_causal <- 4 # number of causal SNPs in region
cov_E_min <- 0
cov_E_max <- 0.5
cov_B_min <- 0.7
cov_B_max <- 0.9

library_size = c(30000000)
read_length = c(75)
paired_end = c("T")

####################################################################################
# load helper functions
####################################################################################

extract_gene_id <- function(transcript_id) {
  # split the string by '|'
  parts <- strsplit(transcript_id, "\\|")[[1]]
  # return the gene ID
  return(parts[2])
}

####################################################################################
# select genes/transcripts to simulate
####################################################################################

set.seed(seed)

if(!dir.exists(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis"))){
  dir.create(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis"),recursive=T)
}

transcripts <- scanFasta("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.transcripts.fa")
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)

transcripts <- transcripts %>%
  mutate(GeneID = sapply(TranscriptID, extract_gene_id))

# remove duplicate transcripts
transcripts <- transcripts %>%
  group_by(TranscriptID) %>%
  filter(Occurrence == 1) %>%
  ungroup()

# read the GTF file
gtf_file <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.annotation.gtf"
gtf <- import(gtf_file)
gtf <- as.data.frame(gtf)

# filter to get only autosomal protein-coding genes
autosomal_genes <- gtf[gtf$seqnames %in% paste0('chr',1:22),]
length(unique(autosomal_genes$gene_id)) # 57623

protein_coding_genes <- autosomal_genes[autosomal_genes$gene_type == "protein_coding", ]
protein_coding_genes <- unique(protein_coding_genes$gene_id)
length(protein_coding_genes) # 19022

# filter transcripts to include only those from protein-coding genes
filtered_transcripts <- transcripts %>%
  filter(GeneID %in% protein_coding_genes)
dim(filtered_transcripts) # 156184

# count the number of isoforms for each gene
isoform_counts <- filtered_transcripts %>%
  group_by(GeneID) %>%
  summarise(IsoformCount = n())
summary(isoform_counts)

isoform_counts <- isoform_counts %>%
  filter(IsoformCount >= min_num_tx)
dim(isoform_counts) # 18983

# get the final list of transcripts for genes with at least four isoforms
final_transcripts <- filtered_transcripts %>%
  filter(GeneID %in% isoform_counts$GeneID)
dim(final_transcripts) # 156184

selected_genes <- protein_coding_genes[protein_coding_genes %in% isoform_counts$GeneID]
length(selected_genes) #18983

if(is.numeric(num_genes_selected)){
  selected_genes <- sample(selected_genes,num_genes_selected)  
}

final_transcripts <- filtered_transcripts %>%
  filter(GeneID %in% selected_genes)
dim(final_transcripts) # 156184      7

# extract start and end positions
# using minimum starting position across all transcripts for start
# using maximum ending position across all transcripts for end
# will define window upstream/downstream of these points for eQTL selection
tx_ranges <- gtf %>%
  filter(gene_id %in% selected_genes) %>%
  group_by(gene_id) %>%
  summarise(start = min(start), end = max(end), seqnames = unique(seqnames), strand = unique(strand))

out_df <- merge(x=final_transcripts, y=tx_ranges,by.x="GeneID",by.y="gene_id")
out_df <- out_df %>%
  mutate(chr = gsub("^chr", "", seqnames))
out_df$chr <- as.numeric(out_df$chr)

out_df <- out_df %>%
  arrange(start) %>%
  arrange(chr)


write.table(out_df, file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/anno_selected_genes.txt"),
  sep="\t",quote=F,row.names=F,col.names=T)

# save gene count by chromosome

# select the first instance of each unique GeneID
unique_genes <- out_df %>%
  group_by(GeneID) %>%
  slice(1) %>%  # grab the first row for each unique GeneID
  ungroup()     

genes_df <- data.frame(count=table(unique_genes$chr))

write.table(genes_df, file=paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/anno_selected_genes_count_by_chr.txt"),
  sep="\t",quote=F,row.names=F,col.names=T)

####################################################################################
# set parameter space
####################################################################################

set.seed(seed)

gene_dat <- data.frame(unique_genes[,c("GeneID","Length","chr","start","end","seqnames","strand")])
gene_dat$prop_shared <- sample(c(0.5, 0.75, 1),nrow(gene_dat),replace=T)
gene_dat$h2_g <- rbeta(nrow(gene_dat), shape1 = 2, shape2 = 40) # changed this here to 40

gene_dat$n_causal <- n_causal
gene_dat$cov_E_min=cov_E_min
gene_dat$cov_E_max=cov_E_max
gene_dat$cov_B_min=cov_B_min
gene_dat$cov_B_max=cov_B_max

# parameter_space <- expand.grid(n_causal = n_causal,
#                               prop_shared = prop_shared,
#                               h2_g = h2_g,
#                               cov_E_min=cov_E_min,
#                               cov_E_max=cov_E_max,
#                               cov_B_min=cov_B_min,
#                               cov_B_max=cov_B_max)

out_dir <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/")
if(!dir.exists(out_dir)){
  dir.create(out_dir,recursive=T)
}

write.table(gene_dat, file=paste0(out_dir,"gene_level_parameters.txt"),
  sep="\t",quote=F,col.names=T,row.names=F)

# set read parameter space

parameter_space <- expand.grid(library_size = library_size,
                              read_length = read_length,
                              paired_end = paired_end)


write.table(parameter_space, file=paste0(out_dir,"parameter_space_reads.txt"),
  sep="\t",quote=F,col.names=T,row.names=F)
