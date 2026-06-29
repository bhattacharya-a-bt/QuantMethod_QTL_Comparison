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
library(data.table)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
pass <- as.numeric(args[1])
param_row_reads <- as.numeric(args[2])

####################################################################################
# begin code
####################################################################################

stats <- data.frame(fread(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/fastqc/param_row_reads_",param_row_reads,"/multiqc_data/multiqc_fastqc.txt")))

# check for failed samples in any column of quality metrics
# filter samples where all relevant columns have "pass"
metrics <- c(
  "basic_statistics", "per_base_sequence_quality", "per_sequence_quality_scores",
  "per_base_sequence_content", "per_sequence_gc_content", "per_base_n_content",
  "sequence_length_distribution", "sequence_duplication_levels",
  "overrepresented_sequences", "adapter_content"
)
cols <- which(names(stats)%in%metrics)
sub <- stats[,cols]
passing <- apply(sub,1,FUN=function(x) sum(x=="fail"))

# Identify samples where all metrics are "pass"
samples_to_keep <- stats$Sample
if (any(passing > 0)) {
    samples_to_keep <- stats$Sample[-which(passing > 0)]
}

write.table(samples_to_keep, paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass,"/files_for_analysis/samples_pass_qc_prreads_",param_row_reads,".txt"),
 row.names = FALSE, col.names = FALSE, quote = FALSE)