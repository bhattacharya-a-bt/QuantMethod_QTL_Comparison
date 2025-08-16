suppressPackageStartupMessages({
  library(Rsubread)
  library(rtracklayer)
  library(data.table)
  library(optparse)
})

# Argument parsing
option_list = list(
  make_option(c("-i", "--index"), type="numeric", help="Index", metavar="numeric"),
  make_option(c("-v", "--version"),type="numeric",help='Version',metavar='numeric')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$index)) {
  stop("Please provide a tissue index with --index")
}

index = opt$index
version = opt$version



SAMPLES_FILE="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/bam_files_with_sample_attrib_rnaseq.txt"
samples = data.table::fread(SAMPLES_FILE)
tissue = unique(samples$SMTSD)[index]

# # build index
# dir.create('/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/rsubread_GRCh38')
setwd('/rsrch5/scratch/epi/bhattacharya_lab/rsubread_GRCh38')
# ref = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta" 
# buildindex(basename="hg38",reference=ref)

if (version == 45){
# assign reads to protein coding genes
## annotation: GeneID Chr Start  End Strand
gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"
gtf.df = data.frame(import(gtf))
annot = gtf.df[gtf.df$type=="gene",]
annot = annot[,c("gene_id","seqnames","start","end","strand")]
names(annot) = c("GeneID","Chr","Start","End","Strand")

DIR_BAM="/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam"
tissue_samples = subset(samples,
                        SMTSD == tissue)
bams = file.path(DIR_BAM,tissue_samples$filename)

if (!file.exists(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45_star',
                           paste0(tissue,'_v45.RDS')))){
counts_v45 = featureCounts(bams,annot.ext=annot,isPairedEnd = TRUE,nthreads=10) # check documentation RE: multi-mapping and multi-feature mapping
saveRDS(list(annot = annot,
             counts = counts_v45),
        file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45_star',
                  paste0(tissue,'_v45.RDS')))}
}

if (version == 38){
# assign reads to protein coding genes
## annotation: GeneID Chr Start  End Strand
gtf = "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.annotation.gtf"
gtf.df = data.frame(import(gtf))
annot = gtf.df[gtf.df$type=="gene",]
annot = annot[,c("gene_id","seqnames","start","end","strand")]
names(annot) = c("GeneID","Chr","Start","End","Strand")
DIR_BAM="/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam"
tissue_samples = subset(samples,
                        SMTSD == tissue)
bams = file.path(DIR_BAM,tissue_samples$filename)

if (!file.exists(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev38_star',
                           paste0(tissue,'_v38.RDS')))){
counts_v38 = featureCounts(bams,annot.ext=annot,isPairedEnd = TRUE,nthreads=10) # check documentation RE: multi-mapping and multi-feature mapping
saveRDS(list(annot = annot,
             counts = counts_v38),
        file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev38_star',
                  paste0(tissue,'_v38.RDS')))}
}

