
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(data.table)
library(cowplot)
library(rtracklayer)

egene_dat <- readRDS("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_eGene_lists.RDS")

egene_dat$setting <- paste(egene_dat$Method,egene_dat$Annotation,egene_dat$Tissue,sep="/")


dat <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_coloc.txt"))
dat <- dat[dat$P<5e-8,]

dat$setting <- paste(dat$quant,dat$annot,dat$tissue,sep="/")
settings <- unique(dat$setting)

names(dat)[names(dat)=="annot"] <- "Annotation"
names(dat)[names(dat)=="quant"] <- "Method"
names(dat)[names(dat)=="tissue"] <- "Tissue"

gwas <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_gwas_lead_snps.txt"))

setDT(dat)
setDT(gwas)

gwas[, CHRPOS := sub("^chr", "", CHRPOS)]

dat[, c("CHR","POS") := tstrsplit(CHRPOS, ":", fixed=TRUE)]
gwas[, c("CHR","POS") := tstrsplit(CHRPOS, ":", fixed=TRUE)]

dat[, POS := as.integer(POS)]
gwas[, POS := as.integer(POS)]

dat[, `:=`(start = POS, end = POS)]
gwas[, `:=`(
  start = POS - 1e6,
  end   = POS + 1e6
)]

setkey(dat, pheno, CHR, start, end)
setkey(gwas, pheno, CHR, start, end)

out <- foverlaps(dat, gwas, nomatch = NA)

setnames(out,
         old = c("ID","POS"),
         new = c("gwas_snp","gwas_pos"))

out <- out[!is.na(out$gwas_snp),]
out$pheno_gwas_snp <- paste0(out$pheno,":",out$gwas_snp)

out$Annotation <- factor(out$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                             labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
out$Method <- factor(out$Method,levels=c("featureCounts","kallisto","salmon","RSEM"),
                         labels=c("featureCounts","kallisto","Salmon","RSEM"))

long_genes <- out %>%
  select(Tissue, Annotation, Method, pheno, gwas_snp) %>%
  unnest(gwas_snp)

long_genes <- unique(long_genes)

# figure 2a_2

gene_overlap<- long_genes %>%
  group_by(Tissue, Method, gwas_snp) %>%
  summarise(
    annotations = list(sort(unique(Annotation))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(annotations, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 annotations")
      if (length(x) == 3) return("3 annotations")
      if (length(x) == 4) return("All 4 annotations")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Tissue, Method, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("GENCODEv27 only","GENCODEv38 only","GENCODEv45 only","Ensembl only","2 annotations","3 annotations","All 4 annotations"))


overlap_prop <- overlap_counts %>%
  group_by(Method, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 annotations",]
# A tibble: 4 × 4
# Groups:   Method [4]
#   Method        overlap_group         n  prop
#   <fct>         <fct>             <int> <dbl>
# 1 featureCounts All 4 annotations 16793 0.425
# 2 kallisto      All 4 annotations 17054 0.420
# 3 Salmon        All 4 annotations 16286 0.408
# 4 RSEM          All 4 annotations 17031 0.446


gg_2a_2 <- ggplot(overlap_prop,
                aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(

  values = c(

    "All 4 annotations" = "#80796b", 

    # single-annotation only groups

    "Ensembl only"     = "#D55E00", 
    "GENCODEv27 only"  = "#E69F00",  
    "GENCODEv38 only"  = "#009E73",  
    "GENCODEv45 only"  = "#CC79A7",  

    # overlap groups

    "2 annotations"    = "#56B4E9",   
    "3 annotations"    = "#0072B2"    

  )

) +
  labs(
    x = "Method", 
    y = "Distribution of GWAS loci\n tagged by colocalization",
    fill = "Annotation"
  )

# figure 2a_1 


gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, gwas_snp) %>%
  summarise(
    methods = list(sort(unique(Method))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 methods")
      if (length(x) == 3) return("3 methods")
      if (length(x) == 4) return("All 4 methods")
    })
  )


overlap_counts <- gene_overlap %>%
  count(Tissue, Annotation, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("featureCounts only","kallisto only","RSEM only","Salmon only","2 methods","3 methods","All 4 methods"))


overlap_prop <- overlap_counts %>%
  group_by(Annotation, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Annotation) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 methods",]
# A tibble: 4 × 4
# Groups:   Annotation [4]
#   Annotation overlap_group     n  prop
#   <fct>      <fct>         <int> <dbl>
# 1 GENCODEv27 All 4 methods 15140 0.373
# 2 GENCODEv38 All 4 methods 15586 0.375
# 3 GENCODEv45 All 4 methods 15466 0.371
# 4 Ensembl    All 4 methods 15870 0.367

gg_2a_1 <- ggplot(overlap_prop,
                aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f",  
"kallisto only"      = "#10b5a2",  
"Salmon only"        = "#b5107b",  
"RSEM only"          = "#6A3D9A",  
"All 4 methods" = "#80796b", 
"2 methods"    = "#56B4E9",   
"3 methods"    = "#0072B2"  
)) +
  labs(
    x = "Annotation", 
    y = "Distribution of GWAS loci\n tagged by colocalization",
    fill = "Method"
  )




# figure 2b_2

dat <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_coloc.txt"))
dat <- dat[dat$P<5e-8,]

names(dat)[names(dat)=="annot"] <- "Annotation"
names(dat)[names(dat)=="quant"] <- "Method"
names(dat)[names(dat)=="tissue"] <- "Tissue"

dat$Annotation <- factor(dat$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                         labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat$Method <- factor(dat$Method,levels=c("featureCounts","kallisto","salmon","RSEM"),
                     labels=c("featureCounts","kallisto","Salmon","RSEM"))

long_genes <- dat %>%
  select(Tissue, Annotation, Method, pheno, phe_id) %>%
  unnest(phe_id)

long_genes <- unique(long_genes)

gene_overlap <- long_genes %>%
  group_by(Tissue, Method, phe_id) %>%
  summarise(
    methods = list(sort(unique(Annotation))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 annotations")
      if (length(x) == 3) return("3 annotations")
      if (length(x) == 4) return("All 4 annotations")
    })
  )  

overlap_counts <- gene_overlap %>%
  count(Tissue, Method, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("GENCODEv27 only","GENCODEv38 only","GENCODEv45 only","Ensembl only","2 annotations","3 annotations","All 4 annotations"))

overlap_prop <- overlap_counts %>%
  group_by(Method, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 annotations",]
# # A tibble: 4 × 4
# # Groups:   Method [4]
#   Method        overlap_group         n  prop
#   <fct>         <fct>             <int> <dbl>
# 1 featureCounts All 4 annotations  9815 0.271
# 2 kallisto      All 4 annotations  9943 0.259
# 3 Salmon        All 4 annotations  9527 0.254
# 4 RSEM          All 4 annotations 10033 0.290

gg_2b_2 <- ggplot(overlap_prop,
                aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("All 4 annotations" = "#80796b",  # dark gray
"Ensembl only"     = "#D55E00",   
"GENCODEv27 only"  = "#E69F00",  
"GENCODEv38 only"  = "#009E73",  
"GENCODEv45 only"  = "#CC79A7",   
"2 annotations"    = "#56B4E9",   
"3 annotations"    = "#0072B2"   
)) +
  labs(
    x = "Method", 
    y = "Proportion of Colocalized Genes",
    fill = "Annotation"
  )

# figure 2b_1

gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, phe_id) %>%
  summarise(
    methods = list(sort(unique(Method))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 methods")
      if (length(x) == 3) return("3 methods")
      if (length(x) == 4) return("All 4 methods")
    })
  )  

overlap_counts <- gene_overlap %>%
  count(Tissue, Annotation, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("featureCounts only","kallisto only","RSEM only","Salmon only","2 methods","3 methods","All 4 methods"))

overlap_prop <- overlap_counts %>%
  group_by(Annotation, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Annotation) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 methods",]
# # A tibble: 4 × 4
# # Groups:   Annotation [4]
#   Annotation overlap_group     n  prop
#   <fct>      <fct>         <int> <dbl>
# 1 GENCODEv27 All 4 methods  8602 0.228
# 2 GENCODEv38 All 4 methods  8876 0.229
# 3 GENCODEv45 All 4 methods  8876 0.225
# 4 Ensembl    All 4 methods  8866 0.210

gg_2b_1 <- ggplot(overlap_prop,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f", 
"kallisto only"      = "#10b5a2",  
"Salmon only"        = "#b5107b",  
"RSEM only"          = "#6A3D9A",   
"All 4 methods" = "#80796b",  
"2 methods"    = "#56B4E9",  
"3 methods"    = "#0072B2"   
))+
  labs(
    x = "Annotation", 
    y = "Proportion of Colocalized Genes",
    fill = "Method"
  )





## figure 2c_2

dat <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_twas_z.txt"))
dat$P <- 2 * pnorm(-abs(dat$TWAS_Z))
dat <- dat[dat$P<2.5e-06,] 

dat$Annotation <- factor(dat$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                             labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat$Method <- factor(dat$Method,levels=c("featureCounts","kallisto","RSEM","salmon"),
                         labels=c("featureCounts","kallisto","Salmon","RSEM"))

gtf_gencode <- data.frame(import("/Users/sthead/OneDrive - Inside MD Anderson/isoqtl_GTEx/pass1/files_for_analysis/annot/gencode/gencode.v45.annotation.gtf"))
gtf_genes <- gtf_gencode[gtf_gencode$type == "gene", ]
gtf_genes$Gene <- sub("\\..*", "", gtf_genes$gene_id)
gtf_genes_small <- gtf_genes[, c("Gene", "seqnames", "start", "end", "strand")]
colnames(gtf_genes_small) <- c("Gene", "CHR", "START", "END", "STRAND")

dat_merged <- merge(dat, gtf_genes_small, by = "Gene", all.x = TRUE,all.y=F)
dat <- dat_merged

gwas <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_gwas_lead_snps.txt"))
gwas <- gwas[gwas$pheno %in% dat$Phenotype,]

setDT(dat)
setDT(gwas)

gwas[, CHRPOS := sub("^chr", "", CHRPOS)]
dat[, CHR := sub("^chr", "", CHR)]
gwas$CHR <- gwas$CHROM

dat[, `:=`(
  start = START - 1e6,
  end   = END + 1e6
)]

dat[start < 0, start := 0]   # prevent negative positions

gwas[, c("CHR","POS") := tstrsplit(CHRPOS, ":", fixed=TRUE)]
gwas[, POS := as.integer(POS)]

gwas[, `:=`(
  start = POS,
  end   = POS
)]

setkey(dat, Phenotype, CHR, start, end)
setkey(gwas, pheno, CHR, start, end)

dat <- dat[!is.na(CHR)]
out <- foverlaps(gwas, dat, nomatch = 0)

colnames(out)[colnames(out)=="ID"] <- "gwas_snp"
long_genes <- out %>%
  select(Tissue, Annotation, Method, pheno, gwas_snp) %>%
  unnest(gwas_snp)

long_genes <- unique(long_genes)

gene_overlap<- long_genes %>%
  group_by(Tissue, Method, gwas_snp) %>%
  summarise(
    annotations = list(sort(unique(Annotation))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(annotations, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 annotations")
      if (length(x) == 3) return("3 annotations")
      if (length(x) == 4) return("All 4 annotations")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Tissue, Method, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("GENCODEv27 only","GENCODEv38 only","GENCODEv45 only","Ensembl only","2 annotations","3 annotations","All 4 annotations"))

overlap_prop <- overlap_counts %>%
  group_by(Method, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 annotations",]
# A tibble: 4 × 4
# Groups:   Method [4]
#   Method        overlap_group         n  prop
#   <fct>         <fct>             <int> <dbl>
# 1 featureCounts All 4 annotations 66677 0.730
# 2 kallisto      All 4 annotations 65907 0.705
# 3 Salmon        All 4 annotations 65475 0.725
# 4 RSEM          All 4 annotations 66215 0.716  

gg_2c_2 <- ggplot(overlap_prop,
                  aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(

  values = c(

    "All 4 annotations" = "#80796b", 

    # single-annotation only groups

    "Ensembl only"     = "#D55E00",  
    "GENCODEv27 only"  = "#E69F00",   
    "GENCODEv38 only"  = "#009E73",  
    "GENCODEv45 only"  = "#CC79A7",   

    # overlap groups

    "2 annotations"    = "#56B4E9",   
    "3 annotations"    = "#0072B2"  

  )

) +
  labs(
    x = "Method", 
    y = "Distribution of GWAS loci\n tagged by TWAS",
    fill = "Annotation"
  )

# figure 2_c_1

gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, gwas_snp) %>%
  summarise(
    methods = list(sort(unique(Method))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 methods")
      if (length(x) == 3) return("3 methods")
      if (length(x) == 4) return("All 4 methods")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Tissue, Annotation, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("featureCounts only","kallisto only","RSEM only","Salmon only","2 methods","3 methods","All 4 methods"))

overlap_prop <- overlap_counts %>%
  group_by(Annotation, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Annotation) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 methods",]
# A tibble: 4 × 4
# Groups:   Annotation [4]
#   Annotation overlap_group     n  prop
#   <fct>      <fct>         <int> <dbl>
# 1 GENCODEv27 All 4 methods 64919 0.699
# 2 GENCODEv38 All 4 methods 65349 0.705
# 3 GENCODEv45 All 4 methods 66314 0.712
# 4 Ensembl    All 4 methods 65249 0.707

gg_2c_1 <- ggplot(overlap_prop,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f",  
"kallisto only"      = "#10b5a2",  
"Salmon only"        = "#b5107b",  
"RSEM only"          = "#6A3D9A",   
"All 4 methods" = "#80796b",  
"2 methods"    = "#56B4E9",   
"3 methods"    = "#0072B2"   
)) +
  labs(
    x = "Annotation", 
    y = "Distribution of GWAS loci\n tagged by TWAS",
    fill = "Method"
  )


# figure 2d_2

dat <- data.frame(fread("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_twas_z.txt"))
dat$P <- 2 * pnorm(-abs(dat$TWAS_Z))
dat <- dat[dat$P<2.5e-06,] 

dat$Annotation <- factor(dat$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                         labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat$Method <- factor(dat$Method,levels=c("featureCounts","kallisto","salmon","RSEM"),
                     labels=c("featureCounts","kallisto","Salmon","RSEM"))


dat$phe_id <- dat$Gene
dat$pheno <- dat$Phenotype

long_genes <- dat %>%
  select(Tissue, Annotation, Method, pheno, phe_id) %>%
  unnest(phe_id)

long_genes <- unique(long_genes)

gene_overlap<- long_genes %>%
  group_by(Tissue, Method, phe_id) %>%
  summarise(
    annotations = list(sort(unique(Annotation))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(annotations, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 annotations")
      if (length(x) == 3) return("3 annotations")
      if (length(x) == 4) return("All 4 annotations")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Tissue, Method, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("GENCODEv27 only","GENCODEv38 only","GENCODEv45 only","Ensembl only","2 annotations","3 annotations","All 4 annotations"))

overlap_prop <- overlap_counts %>%
  group_by(Method, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 annotations",]
# A tibble: 4 × 4
# Groups:   Method [4]
#   Method        overlap_group         n  prop
#   <fct>         <fct>             <int> <dbl>
# 1 featureCounts All 4 annotations 66291 0.247
# 2 kallisto      All 4 annotations 64807 0.219
# 3 Salmon        All 4 annotations 63088 0.224
# 4 RSEM          All 4 annotations 65586 0.259  

gg_2d_2 <- ggplot(overlap_prop,
                  aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("All 4 annotations" = "#80796b",  
"Ensembl only"     = "#D55E00",  
"GENCODEv27 only"  = "#E69F00",  
"GENCODEv38 only"  = "#009E73",
"GENCODEv45 only"  = "#CC79A7",   
"2 annotations"    = "#56B4E9",   
"3 annotations"    = "#0072B2"   
)) +
  labs(
    x = "Method", 
    y = "Proportion of TWAS Genes",
    fill = "Annotation"
  )


# figure 2d_1

gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, phe_id) %>%
  summarise(
    methods = list(sort(unique(Method))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("2 methods")
      if (length(x) == 3) return("3 methods")
      if (length(x) == 4) return("All 4 methods")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Tissue, Annotation, overlap_group)

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,
  levels=c("featureCounts only","kallisto only","RSEM only","Salmon only","2 methods","3 methods","All 4 methods"))

overlap_prop <- overlap_counts %>%
  group_by(Annotation, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Annotation) %>%
  mutate(prop = n / sum(n))

overlap_prop[overlap_prop$overlap_group=="All 4 methods",]
# A tibble: 4 × 4
# Groups:   Annotation [4]
#   Annotation overlap_group     n  prop
#   <fct>      <fct>         <int> <dbl>
# 1 GENCODEv27 All 4 methods 57898 0.208
# 2 GENCODEv38 All 4 methods 59449 0.213
# 3 GENCODEv45 All 4 methods 59492 0.206
# 4 Ensembl    All 4 methods 59740 0.194

gg_2d_1 <- ggplot(overlap_prop,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f",  
"kallisto only"      = "#10b5a2", 
"Salmon only"        = "#b5107b",  
"RSEM only"          = "#6A3D9A",  
"All 4 methods" = "#80796b",  
"2 methods"    = "#56B4E9",   
"3 methods"    = "#0072B2"   
))+
  labs(
    x = "Annotation", 
    y = "Proportion of TWAS Genes",
    fill = "Method"
  )




# cowplot and print to PDF


fin_fig <- cowplot::plot_grid(
  gg_2a_2, gg_2a_1, gg_2b_2, gg_2b_1,
  gg_2c_2, gg_2c_1, gg_2d_2, gg_2d_1,
  labels =c("A","","B","","C","","D"),
  label_x=0,
  label_y=1,
  nrow = 2,
  ncol = 4,
  rel_widths=c(1,1,1,1)
)

ggsave(
  filename = "/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/Manuscript/Figures/r1/r1_fig2_all_annots_and_methods.pdf",
  plot = fin_fig,
  width = 15,     
  height = 9,     
  units = "in"
)