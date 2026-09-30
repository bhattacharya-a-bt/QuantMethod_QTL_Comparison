
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(data.table)
library(cowplot)
library(rtracklayer)

egene_dat <- readRDS("figure_r1_aggregated_eGene_lists.RDS")

egene_dat$setting <- paste(egene_dat$Method,egene_dat$Annotation,egene_dat$Tissue,sep="/")


dat <- data.frame(fread("figure_r1_aggregated_coloc.txt"))
dat <- dat[dat$P<5e-8,]

dat$setting <- paste(dat$quant,dat$annot,dat$tissue,sep="/")
settings <- unique(dat$setting)

names(dat)[names(dat)=="annot"] <- "Annotation"
names(dat)[names(dat)=="quant"] <- "Method"
names(dat)[names(dat)=="tissue"] <- "Tissue"

gwas <- data.frame(fread("figure_r1_gwas_lead_snps.txt"))

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

gg_2a_2_dat <- overlap_prop



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


gg_2a_1_dat <- overlap_prop






# figure 2b_2

dat <- data.frame(fread("figure_r1_aggregated_coloc.txt"))
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


gg_2b_2_dat <- overlap_prop



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


gg_2b_1_dat <- overlap_prop






## figure 2c_2

dat <- data.frame(fread("figure_r1_aggregated_twas_z.txt"))
dat$P <- 2 * pnorm(-abs(dat$TWAS_Z))
dat <- dat[dat$P<2.5e-06,] 

dat$Annotation <- factor(dat$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                             labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat$Method <- factor(dat$Method,levels=c("featureCounts","kallisto","RSEM","salmon"),
                         labels=c("featureCounts","kallisto","Salmon","RSEM"))

gtf_genes_small <- data.frame(fread("figure_gene_info.txt"))

dat_merged <- merge(dat, gtf_genes_small, by = "Gene", all.x = TRUE,all.y=F)
dat <- dat_merged

gwas <- data.frame(fread("figure_r1_gwas_lead_snps.txt"))
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


gg_2c_2_dat <- overlap_prop



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


gg_2c_1_dat <- overlap_prop




# figure 2d_2

dat <- data.frame(fread("figure_r1_aggregated_twas_z.txt"))
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


gg_2d_2_dat <- overlap_prop




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

gg_2d_1_dat <- overlap_prop









#####
########

gg_2a_2 <- ggplot(gg_2a_2_dat,
                aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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

gg_2a_1 <- ggplot(gg_2a_1_dat,
                aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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


gg_2b_2 <- ggplot(gg_2b_2_dat,
                aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
   axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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


gg_2b_1 <- ggplot(gg_2b_1_dat,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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


gg_2c_2 <- ggplot(gg_2c_2_dat,
                  aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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


gg_2c_1 <- ggplot(gg_2c_1_dat,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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

gg_2d_2 <- ggplot(gg_2d_2_dat,
                  aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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

gg_2d_1 <- ggplot(gg_2d_1_dat,
                  aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
axis.text.y  = element_text(size = 6),
axis.title.x = element_text(size = 7),
axis.title.y = element_text(size = 7),
legend.text  = element_text(size = 5),
legend.title = element_text(size = 7),
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


legend_annotation <- get_legend(
  gg_2a_2 +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.box = "horizontal")
)

legend_method <- get_legend(
  gg_2a_1 +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.box = "horizontal")
)

# Hide legends in all panels
plots <- list(
  gg_2a_2, gg_2a_1, gg_2b_2, gg_2b_1,
  gg_2c_2, gg_2c_1, gg_2d_2, gg_2d_1
)
plots <- lapply(plots, function(p) p + theme(legend.position = "none"))

# Rebuild the 2 × 4 plot grid
panel_grid <- plot_grid(
  plotlist = plots,
  labels = c("A", "", "B", "", "C", "", "D", ""),
  label_x = 0, label_y = 1,
  nrow = 2, ncol = 4,
  label_size=7,
  rel_widths = c(1, 1, 1, 1)
)

fin_fig <- plot_grid(
  legend_annotation,
  legend_method,
  panel_grid,
  ncol = 1,
  rel_heights = c(0.04, 0.04, 1)
)


ggsave(
  filename = "figure2.pdf",
  plot = fin_fig,
  width = 180,
  height = 120,
  units = "mm"
)