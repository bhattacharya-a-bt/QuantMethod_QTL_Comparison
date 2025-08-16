#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide a tissue name as the first argument and a numeric index as the second argument.", call. = FALSE)
}
tissue <- args[1]
index <- as.numeric(args[2])

# Check if index conversion was successful
if (is.na(index)) {
  stop("The second argument (index) must be a valid numeric value.", call. = FALSE)
}

cat("Processing tissue:", tissue, "\n")
cat("Index:", index, "\n")
require(data.table)
require(bigsnpr)
require(isotwas)

dir.in = file.path('/rsrch5/scratch/epi/bhattacharya_lab',
                   tissue)

bed_files = list.files(dir.in)
bed_files = bed_files[grepl('bed.gz',bed_files)]
bed_files = bed_files[!grepl('tbi',bed_files)]

for (bbb in bed_files){
  
  setwd(dir.in)
  print(bbb)
  if (!file.exists(file.path(dir.in,
                             paste0(strsplit(bbb,'.bed.gz')[[1]][1],'_withstrand.bed.gz')))){
    
    this = fread(file.path(dir.in,bbb))
    this_output = file.path(dir.in,
                            paste0(strsplit(bbb,'.bed.gz')[[1]][1],'_withstrand.bed'))
    cat('Processing',bbb,'for strand')
    if (!'strnd' %in% colnames(this)){
      ccc = colnames(this)[grepl('GTEX-',colnames(this))]
      this$strnd = '+'
      this = this[,c('#Chr','start','end','pid','gid','strnd',ccc),with=F]
      fwrite(this,
             this_output,
             sep='\t',
             col.names=T,
             row.names=F,
             quote=F)
      
     
      system(paste("bgzip -c", this_output, ">", paste0(this_output, ".gz")))
      system(paste("tabix -p bed", paste0(this_output, ".gz")))}
  }
  
}


temp_dir = file.path('/rsrch5/scratch/epi/abhattacharya3/temp_TWAS',tissue)
dir.create(temp_dir,recursive = T)
VCFFILE="/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/GTEx.WGS.838.vcf.gz"
COVARFILE=paste0("/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/",tissue,"/",tissue,"_formatted_covariates.txt")

bed_files = list.files(dir.in)
bed_files = bed_files[grepl('_withstrand.bed.gz',bed_files)]
bed_files = bed_files[!grepl('tbi',bed_files)]

for (bbb in bed_files){
  cat('Covariate residualizing',bbb)
  infile = file.path(dir.in,bbb)
  outfile = file.path(dir.in,
                      paste0(strsplit(bbb,'.bed.gz')[[1]][1],'_corrected.bed'))
  if (!file.exists(outfile)){
    system(paste('QTLtools correct --bed ',infile,
                 '--cov',COVARFILE,'--normal --out',outfile))}
  
}

bed_files = list.files(dir.in)
names_files <- bed_files[grepl('_corrected.bed',bed_files)]
names = sapply(strsplit(names_files,paste0(tissue,'_')),function(x) x[2])

library(dplyr)
library(stringr)

# Function to parse dataset names
parse_names <- function(names_vector) {
  
  # Create data frame with the names
  df <- data.frame(
    name = names_vector,
    stringsAsFactors = FALSE
  )
  
  # Extract components using regex
  df <- df %>%
    mutate(
      # Extract tissue (everything before the first _v)
      tissue = tissue,
      
      # Extract version (v followed by numbers)
      version = str_extract(name, "v\\d+"),
      
      # Extract method (star or salmon)
      method = str_extract(name, "(star|salmon)"),
      
      # Check if it has "common" in the name
      common = str_detect(name, "common")
    )
  
  return(df)
}


# Parse the names
result_df <- parse_names(names)

# Display the result
print(result_df)


# Read each file and assign to variables with the clean names
for(i in seq_along(names_files)) {
  assign(names[i], data.table::fread(file.path(dir.in,
                                               names_files[i])))
}

gtex_snpfile = '/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF'
snps = snp_attach(snp_readBed2(paste0(gtex_snpfile,'.bed'),
                               backingfile = tempfile()))
dir.create(temp_dir,recursive = T)
snps_matrix = as.matrix(snps$genotypes[])
colnames(snps_matrix) = snps$map$marker.ID
rownames(snps_matrix) = snps$fam$sample.ID

# Function to process TWAS models for a single gene (one row)
process_single_gene <- function(row_index, gene_data, version_name, tissue, snps_matrix, snps) {
  
  out_dir = file.path('/rsrch5/home/epi/bhattacharya_lab/data/TWAS_weights',
                      version_name, tissue)
  dir.create(out_dir, recursive = T, showWarnings = FALSE)
  
  cat('Extracting information for', version_name, 'row', row_index, '...\n')
  gene = gene_data$id[row_index]
  chr = gene_data$`#chr`[row_index]
  start = gene_data$start[row_index]
  end = gene_data$end[row_index]
  
  if (!file.exists(file.path(out_dir, paste0(gene, '_TWAS.RDS')))){
    
    snps_matrix_this = snps_matrix[colnames(gene_data)[7:ncol(gene_data)],
                                   snps$map$marker.ID[snps$map$chromosome == chr &
                                                        snps$map$physical.pos < end + 1e6 &
                                                        snps$map$physical.pos > start - 1e6]]
    
    gene_exp_vector = as.numeric(as.data.frame(gene_data)[row_index, 7:ncol(gene_data)])
    
    cat(paste0('Training models for ', gene, ' at Chr', chr, ' ', start, ':', end, '\n'))
    
    enet = univariate_elasticnet(X = snps_matrix_this,
                                 Y = as.matrix(gene_exp_vector),
                                 Omega = NA,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = NULL,
                                 seed = NULL)
    
    blup = univariate_blup(X = snps_matrix_this,
                           Y = as.matrix(gene_exp_vector),
                           Omega = NA,
                           scale = F,
                           nfolds = 5,
                           verbose = F,
                           par = F,
                           n.cores = NULL,
                           tx_names = NULL,
                           seed = NULL)
    
    susie = univariate_susie(X = snps_matrix_this,
                             Y = as.matrix(gene_exp_vector),
                             Omega = NA,
                             scale = F,
                             alpha = 0.5,
                             nfolds = 5,
                             verbose = F,
                             par = F,
                             n.cores = NULL,
                             tx_names = NULL,
                             seed = NULL)
    
    cat('Finding optimal model\n')
    r2 = unlist(c(enet[[1]]$R2,
                  blup[[1]]$R2,
                  susie[[1]]$R2))
    
    if (max(r2) > .01){
      r2.min = which(r2 == max(r2))
      if (r2.min == 1){
        last.mod = enet
        if (nrow(enet[[1]]$Model) == 0){
          r2.min = which(r2 == median(r2))
        }
      }
      
      if (r2.min == 2){
        last.mod = blup
      }
      
      if (r2.min == 3){
        last.mod = susie
      }
      
      gene.model = last.mod[[1]]$Model
      gene.model$Feature = gene
      gene.model$R2 = last.mod[[1]]$R2
      gene.model$marker.ID = gene.model$SNP
      gene.model = merge(gene.model, snps$map, by='marker.ID')
      colnames(gene.model) = c('marker.ID','SNP',
                               'Weight','Feature','R2','Chromosome',
                               'genetic.dist','Position',
                               'ALT','REF')
      gene.model$Build = 'hg38'
      gene.model = gene.model[,c('Feature','SNP','Chromosome','Position','Build',
                                 'ALT','REF','Weight','R2')]
      
      saveRDS(gene.model, file.path(out_dir, paste0(gene, '_TWAS.RDS')))
    } else {
      saveRDS(paste0('R2 < 0.01'), file.path(out_dir, paste0(gene, '_TWAS.RDS')))
    }
  }
  
  return(gene) # Return gene name for tracking
}

# Function to process all rows for a given dataset using lapply
process_twas_models <- function(gene_data, version_name, tissue, snps_matrix, snps) {
  cat("Processing", nrow(gene_data), "genes for", version_name, "...\n")
  
  # Use lapply to process each row with error handling
  results <- lapply(1:nrow(gene_data), function(i) {
    tryCatch({
      process_single_gene(i, gene_data, version_name, tissue, snps_matrix, snps)
    }, error = function(e) {
      cat("Error processing gene at row", i, "in", version_name, ":", e$message, "\n")
      return(paste0("ERROR_ROW_", i))
    })
  })
  
  # Count successful vs failed processes
  errors <- sum(grepl("^ERROR_ROW_", results))
  successes <- length(results) - errors
  
  cat("Completed processing", version_name, ":", successes, "successful,", errors, "errors\n")
  return(results)
}

# Process all datasets from result_df
i = index
  # Get dataset info
  dataset_name <- result_df$name[i]
  tissue <- result_df$tissue[i]
  version <- gsub("v", "", result_df$version[i])  # Remove 'v' prefix
  method <- result_df$method[i]
  common_id <- ifelse(result_df$common[i], "common", "total")
  if (common_id == 'common'){
  # Create version_name for output directory
  version_name <- paste0(version, "_", method, "_", common_id)
  
  # Get the actual data frame using the name
  gene_data <- get(dataset_name)
  
  cat("Processing dataset:", dataset_name, "->", version_name, "\n")
  
  # Process TWAS models
  process_twas_models(gene_data, version_name, tissue, snps_matrix, snps)
  }
