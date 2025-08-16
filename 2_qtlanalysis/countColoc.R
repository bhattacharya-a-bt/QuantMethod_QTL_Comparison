require(data.table)
require(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript countColoc.R <tissue_index> <cancer_index>")
}

tissue_index <- as.numeric(args[1])
cancer_index <- as.numeric(args[2])

setwd('/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare')
tissues = list.files()

# Check if tissue_index is valid
if (tissue_index < 1 || tissue_index > length(tissues)) {
  stop(paste("Invalid tissue_index:", tissue_index, "Valid range: 1 to", length(tissues)))
}

# Get the specific tissue
t <- tissues[tissue_index]
print(paste("Processing tissue:", t, "(index:", tissue_index, ")"))

setwd(file.path(t,'coloc'))
cancers <- list.files()

# Check if cancer_index is valid
if (cancer_index < 1 || cancer_index > length(cancers)) {
  stop(paste("Invalid cancer_index:", cancer_index, "Valid range: 1 to", length(cancers)))
}

# Get the specific cancer
cancer <- cancers[cancer_index]
print(paste("Processing cancer:", cancer, "(index:", cancer_index, ")"))

setwd(cancer)
fff = list.files()

folder = getwd()
if (file.exists(file.path(getwd(),'all_coloc.txt'))){
  file.remove(file.path(getwd(),'all_coloc.txt'))
}
output_file = file.path(getwd(), 'all_coloc.txt')

cat('Aggregating all coloc files\n')
cmd <- sprintf(
  'for f in %s/*.txt.gz; do gzip -cd "$f" | awk -v fname="$(basename "$f")" \'BEGIN{OFS="\\t"} {print $0, fname}\' ; done > %s',
  folder,
  output_file
)
system(cmd)

cat('Processing CLPP\n')
res = fread(output_file)
res = subset(res, phe_id != 'phe_id')
res$CLPP = as.numeric(res$CLPP)
res$P = as.numeric(res$P)
colnames(res)[24] = 'File'
library(dplyr)
library(stringr)

result <- res %>%
  group_by(File) %>%
  summarise(
    Gene = first(phe_id),  # assuming Gene is consistent within each File
    max_CLPP = max(CLPP, na.rm = TRUE),
    min_P = min(P, na.rm = TRUE),
    var_id_min_P = var_id[which.min(P)],
    CHRPOS_min_P = CHRPOS[which.min(P)],
    var_id_max_CLPP = var_id[which.max(CLPP)],
    CHRPOS_max_CLPP = CHRPOS[which.max(CLPP)],
    .groups = "drop"
  ) %>%
  mutate(
    m = map_chr(File, ~ str_split(str_split(.x, "coloc_")[[1]][2], "\\.txt")[[1]][1]),
    Version = map_chr(m, ~ str_split(.x, "_")[[1]][1]),
    Method = map_chr(m, ~ str_split(.x, "_")[[1]][2])
  ) %>%
  select(-m)  # remove the intermediate column if not needed

result$Tissue = t
result$Cancer = cancer
colnames(result) = c('File','Gene','CLPP','P','GWAS_SNP',
                     'GWAS_CHRPOS','CLPP_SNP','CLPP_CHRPOS',
                     'Version','Method','Tissue','Cancer')

result = subset(result, CLPP >= 0.01 & P <= 5e-7)

if (nrow(result) > 0) {
  cat('Writing out!\n')
  
  # Create unique output filename for this tissue-cancer combination
  output_filename <- paste0('colocalization_results_', 
                            gsub('[^A-Za-z0-9_]', '_', t), '_',
                            gsub('[^A-Za-z0-9_]', '_', cancer), '_',
                            tissue_index, '_', cancer_index, '.tsv')
  
  fwrite(result,
         file.path('/rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis', output_filename),
         row.names = FALSE,
         quote = FALSE,
         sep = '\t')
  
  cat(paste('Results written to:', output_filename, '\n'))
} else {
  cat('No hits\n')
}

cat(paste('Completed processing tissue', tissue_index, 'cancer', cancer_index, '\n'))