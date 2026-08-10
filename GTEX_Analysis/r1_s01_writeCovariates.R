####################################################################################
# USER CONFIGURATION - edit these paths, or set the matching environment variables
####################################################################################
# tissue_list_dir  directory whose subdirectory names are the tissues to loop over
# covar_dir        directory holding GTEx '<tissue>.v8.covariates.txt.covar' files
# analysis_dir     output directory; '<tissue>/<tissue>_formatted_covariates.txt'
#                  is written underneath it (created if missing)
tissue_list_dir <- Sys.getenv("TISSUE_LIST_DIR", "/path/to/GTEx_v8/gencodev45")
covar_dir       <- Sys.getenv("COVAR_DIR",       "/path/to/GTEx_v8/covariate_files")
analysis_dir    <- Sys.getenv("ANALYSIS_DIR",    "/path/to/scratch/GTEx_gencode_comp/requant_analyses")
####################################################################################

for (tissue in list.files(tissue_list_dir)){
covar_file <- file.path(covar_dir,
                        paste0(tissue,'.v8.covariates.txt.covar'))
dir.create(file.path(analysis_dir,
                     tissue),
           recursive = T)
out_file <- file.path(analysis_dir,
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