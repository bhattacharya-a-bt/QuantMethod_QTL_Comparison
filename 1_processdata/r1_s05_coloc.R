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
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(bigsnpr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
annot <- as.character(args[1])
quant <- as.character(args[2])
tissue <- as.character(args[3])
pheno_file <- as.character(args[4])
pheno_name <- as.character(args[5])
liftoverTo38 <- as.logical(args[6])
n_bins <- as.integer(args[7])
bin_index <- as.integer(args[8])

####################################################################################
# helper functions
####################################################################################

# function to perform liftover using R liftOver package
liftover_gwas <- function(gwas_data, chain_file) {
  require(rtracklayer)
  require(GenomicRanges)
  require(data.table)
  
  # load chain file
  chain <- import.chain(chain_file)
  
  # create GRanges object
  gwas_gr <- GRanges(
    seqnames = paste0("chr", gwas_data$CHR),
    ranges = IRanges(start = gwas_data$POS, end = gwas_data$POS),
    strand = "*"
  )
  
  # add metadata to track original rows
  mcols(gwas_gr) <- data.frame(
    orig_idx = 1:nrow(gwas_data),
    SNP = gwas_data$SNP,
    CHR = gwas_data$CHR,
    POS = gwas_data$POS
  )
  
  # liftover
  lifted_gr <- liftOver(gwas_gr, chain)
  
  lifted_gr_unlisted <- unlist(lifted_gr)
  
  if (length(lifted_gr_unlisted) > 0) {
    # get original indices of lifted variants
    orig_indices <- mcols(lifted_gr_unlisted)$orig_idx
    
    gwas_lifted <- gwas_data[orig_indices, ]
    
    # update coordinates with lifted positions
    gwas_lifted$CHR <- as.numeric(gsub("chr", "", as.character(seqnames(lifted_gr_unlisted))))
    gwas_lifted$POS <- start(lifted_gr_unlisted)
    
    cat("    Liftover successful:", length(lifted_gr_unlisted), "out of", nrow(gwas_data), "variants lifted\n")
    
    return(gwas_lifted)
  } else {
    warning("Liftover failed for all variants, using original coordinates")
    return(gwas_data)
  }
}

# function to perform ecaviar colocalization analysis
runColoc = function(qtl_data,
                        gwas_name,ldfile,gwas_zfile,
                        gwas_data_lifted){
      
      if (nrow(qtl_data) > 0) {
        qtl_data$CHRPOS <- paste(chr, qtl_data$var_from, sep = ':')
        qtl_data <- subset(qtl_data, CHRPOS %in% int_chrpos)
        
        if (nrow(qtl_data) > 0) {
          zfile <- file.path(tempfolder, paste0('eQTL_',gwas_name,
                                                '.txt'))
          qtl_out <- data.frame(V1 = qtl_data$CHRPOS,
                                V2 = abs(qnorm(qtl_data$nom_pval)) * sign(qtl_data$slope))
          data.table::fwrite(qtl_out, zfile, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          out_ecav <- file.path(tempfolder, paste0('eQTL_',
                                                   gwas_name, '_ecav'))
          system(paste(ecaviar,
                       '-o', out_ecav,
                       '-l', ld_file,
                       '-l', ld_file,
                       '-z', gwas_zfile,
                       '-z', zfile,
                       '-c', 2))
          
          if (file.exists(paste0(out_ecav, '_col'))) {
            ecav <- data.table::fread(paste0(out_ecav, '_col'))
            colnames(ecav)[1] <- 'CHRPOS'
            qtl_data <- merge(qtl_data, ecav, by = 'CHRPOS')
            qtl_data <- merge(qtl_data, gwas_data_lifted[, c('CHRPOS', 'P', 'BETA', 'SE', 'Z')], by = 'CHRPOS')
            
            dir_out <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/coloc_results/",tissue)
            dir.create(dir_out, recursive = TRUE)
            
            data.table::fwrite(qtl_data,
                               file.path(dir_out, paste0(gene, '_eQTL_coloc_',gwas_name,'_',annot,'_',quant,'.txt.gz')),
                               sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
            cat('Written out to',
                file.path(dir_out, paste0(gene, '_eQTL_coloc_',gwas_name,'_',annot,'_',quant,'.txt.gz')))
          }
        }
      }
      
    }

####################################################################################
# begin code
####################################################################################

# gene name and location information
gene_info1 <- data.frame(fread("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/ensembl_gene_info.txt"))

# egenes information
gene_info <- readRDS('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/cis_eqtl_results/r1_aggregated_eGene_lists.RDS')

idx <- which(gene_info$Tissue==tissue & gene_info$Annotation==annot & gene_info$Method==quant)
gene_info <- gene_info[idx,]

egenes <- gene_info$ensembl_gene_id[[1]]
gene_count <- 1:length(egenes)
bins <- cut(
  gene_count,
  breaks = n_bins,
  labels = FALSE
)
egenes = egenes[bins==bin_index]
egenes_info <- gene_info1[gene_info1$ensembl_gene_id %in% egenes,]

# load gwas data
cat(paste0("Reading in GWAS ",pheno_name," \n"))
gwas_data <- data.table::fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/",pheno_file))

# check if columns match expected format and rename if necessary
if(pheno_name %in% c("BreastCancer","ProstateCancer")){
  if ("SNP" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
  # Rename columns to match expected format
  colnames(gwas_data)[colnames(gwas_data) == "BP"] <- "POS"
  colnames(gwas_data)[colnames(gwas_data) == "A1"] <- "ALT"
  colnames(gwas_data)[colnames(gwas_data) == "A2"] <- "REF"
  } else {
    warning("GWAS file format may not match expected columns")
  }
  }else{

if ("SNP" %in% colnames(gwas_data) && "BP" %in% colnames(gwas_data)) {
    # Rename columns to match expected format
    colnames(gwas_data)[colnames(gwas_data) == "BP"] <- "POS"
    colnames(gwas_data)[colnames(gwas_data) == "A1"] <- "REF"
    colnames(gwas_data)[colnames(gwas_data) == "A2"] <- "ALT"
  } else {
    warning("GWAS file format may not match expected columns")
  }
}

# remove prefix chr if present
if(length(grep("^chr",gwas_data$CHR[1]))==1){
    gwas_data[, CHR := sub("^chr", "", CHR)]
    gwas_data$CHR <- as.numeric(gwas_data$CHR)
}

# chain file for liftover
chain_file <- '/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain'

# loop through each gene
for (i in 1:nrow(egenes_info)) {
  
  gene_interest <- egenes_info[i,]
  
  gene <- gene_interest$ensembl_gene_id[1]
  start <- gene_interest$start_position[1]
  chr <- gene_interest$chromosome_name[1]
  end <- gene_interest$end_position[1]
  biotype <- gene_interest$gene_biotype[1]
  
  cat("Processing gene:", gene, "on chromosome", chr, "\n")
  
  gwas_file <- pheno_file
  gwas_name <- pheno_name
  
  cat("  Processing GWAS:", gwas_name, "\n")
  
  
  # Create temp folder for this gene-GWAS combination
  tempfolder <- file.path('/rsrch5/home/epi/sthead/',
                          'tempColoc',
                          paste0('gene', i, '_', gene, '_', gwas_name, "_",tissue,"_",annot,"_",quant))
  dir.create(tempfolder, recursive = TRUE)
  
  
  
  # filter to chromosome region (pre-liftover)
  gwas_data_subset <- subset(gwas_data,
                             CHR == chr)
  
  if (nrow(gwas_data_subset) == 0) {
    cat("    No SNPs found in region for", gwas_name, "\n")
    unlink(tempfolder, recursive = TRUE)
    next
  }
  
  # perform liftover if needed
  if(liftoverTo38==T){
    gwas_data_lifted <- liftover_gwas(gwas_data_subset,chain_file)
  }else{
    gwas_data_lifted <- gwas_data_subset   
  }
   
  # calculate Z-score if not present
  if (!"Z" %in% colnames(gwas_data_lifted)) {
    gwas_data_lifted$Z <- gwas_data_lifted$BETA / gwas_data_lifted$SE
  }
  
  # filter again after liftover to hg38 
  gwas_data_lifted <- subset(gwas_data_lifted,
                             CHR == chr &
                               POS < end + 1e6 &
                               POS > start - 1e6)
  
  if (nrow(gwas_data_lifted) == 0 | min(gwas_data_lifted$P) > 1e-7) { # found a bug here, previously was <= 1e-7
    cat("    No SNPs found in region after liftover for", gwas_name, "\n")
    unlink(tempfolder, recursive = TRUE)
    next
  }
  
  gwas_data_lifted$CHRPOS <- paste(gwas_data_lifted$CHR, gwas_data_lifted$POS, sep = ':')
  
  
  # extract LD from GTEx
  ecaviar <- '/rsrch5/home/epi/bhattacharya_lab/software/caviar/CAVIAR-C++/eCAVIAR'
  gtex_snpfile <- '/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated'
  
  require(bigsnpr)
  system(paste('plink2 --bfile', paste0(gtex_snpfile),
               '--chr', chr,
               '--from-bp', start - 1e6,
               '--to-bp', end + 1e6,
               '--make-bed',
               '--out', file.path(tempfolder, paste0(gene, '_', gwas_name))))
  
  bed_file_path <- file.path(tempfolder, paste0(gene, '_', gwas_name, '.bed'))
  if (!file.exists(bed_file_path)) {
    cat("    Failed to create GTEx bed file for", gwas_name, "\n")
    unlink(tempfolder, recursive = TRUE)
    next
  }
  
  snps_gtex <- snp_attach(snp_readBed2(bed_file_path))
  snps_gtex$map$CHRPOS <- paste(snps_gtex$map$chromosome,
                                snps_gtex$map$physical.pos,
                                sep = ':')
  
  # intersect GWAS AND GTEx
  int_chrpos <- intersect(gwas_data_lifted$CHRPOS, snps_gtex$map$CHRPOS)
  
  if (length(int_chrpos) == 0) {
    cat("    No overlapping SNPs between GWAS and GTEx for", gwas_name, "\n")
    unlink(tempfolder, recursive = TRUE)
    next
  }
  
  gwas_data_lifted <- subset(gwas_data_lifted, CHRPOS %in% int_chrpos)
  snps_gtex <- snp_attach(subset(snps_gtex,
                                 ind.col = which(snps_gtex$map$CHRPOS %in% int_chrpos)))
  gwas_data_lifted <- gwas_data_lifted[match(snps_gtex$map$CHRPOS,
                                             gwas_data_lifted$CHRPOS), ]
  gwas_data_lifted$Z <- ifelse(gwas_data_lifted$REF == snps_gtex$map$allele2,
                               gwas_data_lifted$Z,
                               -1 * gwas_data_lifted$Z)
  
  # calculate LD matrix
  ld <- as.matrix(snp_cor(snps_gtex$genotypes))
  colnames(ld) <- rownames(ld) <- snps_gtex$map$CHRPOS
  ld_file <- file.path(tempfolder, paste0('LD_', gwas_name, '.txt'))
  write.table(ld, ld_file, sep = '\t', col.names = FALSE, row.names = FALSE)
  
  # write GWAS zscores
  gwas_zfile <- file.path(tempfolder, paste0('gwas_', gwas_name, '.txt'))
  gwas_out <- data.frame(V1 = gwas_data_lifted$CHRPOS, V2 = gwas_data_lifted$Z)
  data.table::fwrite(gwas_out, gwas_zfile, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # run qtltools nominal pass

  # first subset bed file to gene of interest
  bed <- (fread(paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/",tissue,"/",annot,"_",quant,".v8.normalized_expression.bed.gz")))
  idxx <- which(bed$pid==gene)
  bed_this <- bed[idxx,]
  write.table(bed_this,file=paste0(tempfolder,"/tmp.bed"),col.names=T,row.names=F,sep="\t",quote=F)

  # create subset vcf file for qtltools
  prefix   <- paste0(gene, "_", gwas_name)
  vcf_file <- paste0(prefix, ".vcf.gz")

  cmd <- paste(
    "cd", shQuote(tempfolder), "&&",
    "plink2",
    "--bfile", shQuote(prefix),
    "--export vcf bgz",
    "--output-chr chrM",
    "--out", shQuote(prefix), "&&",
    "tabix -p vcf", shQuote(vcf_file)
  )

  system(cmd)

  # run qtltools nominal pass
  cov_file <- sprintf(
    "/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/%s/%s_formatted_covariates.txt",
    tissue, tissue
  )

  geno_file=paste0(gene,"_",gwas_name,".vcf.gz")

  cmd <- sprintf(
    "cd %s && bgzip -f tmp.bed && tabix -p bed tmp.bed.gz && QTLtools cis --vcf %s --bed tmp.bed.gz --nominal 1 --std-err --cov %s --out nominal_qtltools_res.txt --normal --seed 1219",
    shQuote(tempfolder),
    shQuote(geno_file),
    shQuote(cov_file)
  )

  system(cmd)

  # moving on to running ecaviar

  header_names <- read.table('/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/nominal_header_0.txt',
                             header = FALSE)
  header_names <- as.character(header_names[1, ])
  header_names = c(header_names,'best_hit')
  header_names[15] = 'std_err'

  #setwd(tempfolder)
  
  # read in qtl data
  qtl_data <- fread(paste0(tempfolder,"/nominal_qtltools_res.txt"))
  colnames(qtl_data) <- header_names

  # run colocalization 
  cat("Beginning colocalization with eCAVIAR... \n")

  runColoc(qtl_data = qtl_data,
           gwas_name = gwas_name,
           ldfile = ld_file,
           gwas_zfile = gwas_zfile,
           gwas_data_lifted = gwas_data_lifted)
      
  # Clean up temp folder for this GWAS
  unlink(tempfolder, recursive = TRUE)
  if (!dir.exists(tempfolder)) {
    cat("Temp folder", tempfolder, "was deleted!\n")
  }
    
  cat("Completed gene:", gene, "\n\n")
  
} # end gene loop

