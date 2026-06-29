for (tissue in list.files('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45')){
covar_file <- file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/covariate_files',
                        paste0(tissue,'.v8.covariates.txt.covar'))
dir.create(file.path('/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/',
                     tissue),
           recursive = T)
out_file <- file.path('/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/',
                      tissue,
                      paste0(tissue,"_formatted_covariates.txt"))

cat("Formatting covariates from: ", covar_file, "\n")

library(dplyr)

# read in data
covariates <- readr::read_delim(file = covar_file,
                                col_names=FALSE)[,-1] # drop redundant first column

# reformat standard GTEx covariates per tissue
ids <- covariates %>%
  pull(X2) #%>%
  #stringr::str_replace_all("\\.", "-")
covariates <- covariates[,-1] %>% t
colnames(covariates) <- ids
covariates <- cbind("id" = c(paste0("V", 1:nrow(covariates))),
                    covariates)
rownames(covariates) <- NULL

# write tab-delimited formatted_covariates.txt
write.table(covariates,
            file = out_file,
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
}