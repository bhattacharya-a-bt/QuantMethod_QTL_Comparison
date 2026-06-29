pass = 2
n = 500

####################################################################################

directory <- paste0("/rsrch5/scratch/epi/sthead/isoqtl_quant_error/pass",pass,"/files_for_analysis")

# check if the directory doesn't exist
if (!file.exists(directory)) {
  # Create the directory along with any intermediate directories
  dir.create(directory, recursive = TRUE)
  cat("Directory created:", directory, "\n")
} else {
  cat("Directory already exists:", directory, "\n")
}

geno_dir = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG"
proj_dir = paste0("/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass",pass)
setwd(proj_dir)

dat <- read.csv(paste0(geno_dir,"/all_hg38.psam"),sep="\t")
eur <- dat[dat$SuperPop=="EUR",]
table(eur$Population)
dim(eur) # 633x6

# read in 1KG metadata from "kgp" R package
# cannot install this package on HPC due to R version
load("files_for_analysis/kgp3.RData")
dim(kgp3) # 2504 samples sequenced for the phase 3 release

eur <- eur[eur$X.IID %in% kgp3$id,]
dim(eur) # 503x6

# randomly select 500
set.seed(12345)
sel <- sample(1:nrow(eur),n)
samp <- data.frame(ID=eur[sel,1])

write.table(samp,file="files_for_analysis/1kg_eur_500_sample_ids",
            sep="\t",row.names = F,col.names = F,quote=F)