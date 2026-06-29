#!/usr/bin/env Rscript

##################################################################################
# change library to local
##################################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

####################################################################################
# load dependencies
####################################################################################
require(edgeR)
require(dplyr)
library(data.table)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
annot <- as.character(args[1])
quant <- as.character(args[2])
tissue <- as.character(args[3])

# save ensembl gene information (run locally once before generating bed files)
# library(biomaRt)
# mart <- useEnsembl(
#   biomart = "genes",
#   dataset = "hsapiens_gene_ensembl"
# )

# # query gene coordinates
# all_genes <- getBM(
#   attributes = c(
#     "ensembl_gene_id",
#     "chromosome_name",
#     "start_position",
#     "end_position",
#     "strand",
#     "gene_biotype",
#     "hgnc_symbol"
#   ),
#   mart = mart
# )
# write.table(all_genes,file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/ensembl_gene_info.txt",
#   sep="\t",quote=F,row.names=F,col.names=T)

####################################################################################
# begin code
####################################################################################

gene_info <- data.frame(fread("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/ensembl_gene_info.txt"))

setwd(paste0('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/polished/',annot))
t = paste0(tissue,"_",quant,"_gene.RDS")

out_file <- paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/",tissue,"/",annot,"_",quant,".v8.normalized_expression.bed")

print(t)
all_in = readRDS(t)

# get CPM with normalization factors applied
tmm_normalized = as.data.frame(all_in$TMM_normalized)

# apply log2 transformation to normalized data
expression_data <- log2(tmm_normalized + 1)

rrr = rownames(tmm_normalized)[which(rowMeans(tmm_normalized <= .1) <= .25)]
expression_data <- expression_data[rrr,]

# some datasets have gene version numbers so need to get rid conditionally

obs_gene <- data.frame(
  "ensembl_gene_id" = sub("\\..*$", "", rownames(expression_data))
)

#obs_gene <- data.frame("ensembl_gene_id"=rownames(expression_data))
bed_out <- cbind(obs_gene,expression_data)

bed_out <- merge(bed_out,gene_info[,c("ensembl_gene_id","chromosome_name","start_position","end_position","strand")],
  by="ensembl_gene_id",all.x=T,all.y=F)

bed_out <- bed_out[,c("chromosome_name","start_position","end_position","ensembl_gene_id","ensembl_gene_id","strand",colnames(bed_out)[grep("^GTEX",colnames(bed_out))])]

colnames(bed_out)[1:6] = c('#Chr','start','end','pid','gid','strand')
bed_out$strand <- as.character(bed_out$strand)
bed_out$strand[bed_out$strand=="1"] <- "+"
bed_out$strand[bed_out$strand=="-1"] <- "-"

bed_out <- data.table(bed_out)

bed_out_sorted <- bed_out %>%
  mutate(
    chr_ord = case_when(
      `#Chr` %in% as.character(1:22) ~ as.integer(`#Chr`),
      `#Chr` == "X" ~ 23L,
      `#Chr` == "Y" ~ 24L,
      TRUE ~ 99L
    )
  ) %>%
  arrange(chr_ord, start) %>%
  select(-chr_ord)



newcols <- sub(
"^(GTEX-[^-]+).*$",
"\\1",
colnames(bed_out_sorted)
)
colnames(bed_out_sorted) <- newcols

colnames(bed_out_sorted) <- gsub("-", ".", colnames(bed_out_sorted))

bed_out_sored <- data.frame(bed_out_sorted)
bed_out_sorted$`#Chr` <- paste0("chr",bed_out_sorted$`#Chr`)

data.table::fwrite(bed_out_sorted,
                   paste0(out_file),
                   sep='\t',
                   col.names=T,
                   row.names=F)



