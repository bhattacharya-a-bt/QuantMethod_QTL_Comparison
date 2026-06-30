
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(data.table)
library(cowplot)

# load data for for top row 

dat_sub <- readRDS("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_eGene_lists.RDS")

dat_sub$Annotation <- factor(dat_sub$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45","Ensembl"),
                             labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat_sub$Method <- factor(dat_sub$Method,levels=c("featureCounts","kallisto","RSEM","salmon"),
                         labels=c("featureCounts","kallisto","Salmon","RSEM"))

long_genes <- dat_sub %>%
  select(Tissue, Annotation, Method, ensembl_gene_id) %>%
  unnest(ensembl_gene_id)

# figure 1b

gene_overlap<- long_genes %>%
  group_by(Tissue, Method, ensembl_gene_id) %>%
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
  
overlap_prop[grep("All 4",overlap_prop$overlap_group),]
# # A tibble: 4 × 4
# # Groups:   Method [4]
#   Method        overlap_group          n  prop
#   <fct>         <chr>              <int> <dbl>
# 1 featureCounts All 4 annotations 284303 0.655
# 2 kallisto      All 4 annotations 287644 0.621
# 3 Salmon        All 4 annotations 280162 0.680
# 4 RSEM          All 4 annotations 274520 0.615

gg_1b <- ggplot(overlap_prop,
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
    y = "Proportion of eGenes",
    fill = "Annotation"
  )

# figure 1d

N_GTEx <- read.delim("~/Library/CloudStorage/OneDrive-InsideMDAnderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/N_GTEx.tsv")

N_GTEx$Tissue[N_GTEx$Tissue=="Brain_Spinal_cord_cervical_c-1"] <- "Brain_Spinal_cord_cervical_c_1"
N_GTEx$Tissue[N_GTEx$Tissue=="Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV_transformed_lymphocytes"

method_consistency <- overlap_counts %>%
  group_by(Tissue, Method) %>%
  summarise(
    n_all4 = sum(n[overlap_group == "All 4 annotations"]),
    total_egenes = sum(n),
    prop_all4 = n_all4 / total_egenes,
    .groups = "drop"
  )

method_consistency <- method_consistency %>%
  mutate(
    discordance = 1 - prop_all4
  )

plotdat <- merge(x=method_consistency, y=N_GTEx,by="Tissue",all.x = T,all.y=F)


gg_1d <- ggplot(plotdat,
                aes(x = N, y = discordance, color = Method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid")  +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),
    legend.position = c(0.02, 1),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "eGene discordance %\nacross annotation",
    color = "Method"
  ) +
  scale_y_continuous(labels = scales::percent_format())+
  scale_color_manual(
  values = c(
    "featureCounts" = "#b9d47f",  
    "kallisto"      = "#10b5a2",  
    "Salmon"        = "#b5107b",  
    "RSEM"          = "#6A3D9A"   
  )
)

# figure 1a

gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, ensembl_gene_id) %>%
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

overlap_prop[grep("All 4",overlap_prop$overlap_group),]
#   Annotation overlap_group      n  prop
#   <fct>      <fct>          <int> <dbl>
# 1 GENCODEv27 All 4 methods 255833 0.597
# 2 GENCODEv38 All 4 methods 266227 0.603
# 3 GENCODEv45 All 4 methods 267288 0.589
# 4 Ensembl    All 4 methods 279947 0.567

gg_1a <- ggplot(overlap_prop,
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
    y = "Proportion of eGenes",
    fill = "Method"
  )

# figure 1c

N_GTEx <- read.delim("~/Library/CloudStorage/OneDrive-InsideMDAnderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/N_GTEx.tsv")

N_GTEx$Tissue[N_GTEx$Tissue=="Brain_Spinal_cord_cervical_c-1"] <- "Brain_Spinal_cord_cervical_c_1"
N_GTEx$Tissue[N_GTEx$Tissue=="Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV_transformed_lymphocytes"

method_consistency <- overlap_counts %>%
  group_by(Tissue, Annotation) %>%
  summarise(
    n_all4 = sum(n[overlap_group == "All 4 methods"]),
    total_egenes = sum(n),
    prop_all4 = n_all4 / total_egenes,
    .groups = "drop"
  )

method_consistency <- method_consistency %>%
  mutate(
    discordance = 1 - prop_all4
  )

plotdat <- merge(x=method_consistency, y=N_GTEx,by="Tissue",all.x = T,all.y=F)


gg_1c <- ggplot(plotdat,
                aes(x = N, y = discordance, color = Annotation)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid")  +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),
    legend.position = c(0.02, 1),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "eGene discordance %\nacross method",
    color = "Annotation"
  ) +
  scale_y_continuous(labels = scales::percent_format())+
  scale_color_manual(

  values = c(

    "Ensembl"     = "#D55E00",   
"GENCODEv27"  = "#E69F00",   
"GENCODEv38"  = "#009E73",   
"GENCODEv45"  = "#CC79A7"  
    ))



# bottom row (predicted TWAS genes)


dat <- readRDS("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_TWAS_passR2.RDS")
colnames(dat)[colnames(dat)=="tissue"] <- "Tissue"
colnames(dat)[colnames(dat)=="annot"] <- "Annotation"
colnames(dat)[colnames(dat)=="quant"] <- "Method"
colnames(dat)[colnames(dat)=="pass_genes"] <- "ensembl_gene_id"

dat_sub <- dat %>%
  mutate( 
    ensembl_gene_id = map(ensembl_gene_id, 1)
  )

dat_sub$Annotation <- factor(dat_sub$Annotation,levels=c("GENCODE_V27","GENCODE_V38","GENCODE_V45","Ensembl"),
                             labels=c("GENCODEv27","GENCODEv38","GENCODEv45","Ensembl"))
dat_sub$Method <- factor(dat_sub$Method,levels=c("featureCounts","kallisto","RSEM","salmon"),
                         labels=c("featureCounts","kallisto","RSEM","Salmon"))

long_genes <- dat_sub %>%
  mutate(ensembl_gene_id = as.list(ensembl_gene_id)) %>% 
  unnest(ensembl_gene_id)

# figure 1b

gene_overlap<- long_genes %>%
  group_by(Tissue, Method, ensembl_gene_id) %>%
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
#   Method        overlap_group          n  prop
#   <fct>         <fct>              <int> <dbl>
# 1 featureCounts All 4 annotations 343919 0.375
# 2 kallisto      All 4 annotations 347180 0.329
# 3 RSEM          All 4 annotations 331230 0.388
# 4 Salmon        All 4 annotations 334422 0.345

gg_1b_bottom <- ggplot(overlap_prop,
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
))+
  labs(
    x = "Method", 
    y = expression("Proportion of genes with"~R^2 >= 0.01),
    fill = "Annotation"
  )

# figure 1d

N_GTEx <- read.delim("~/Library/CloudStorage/OneDrive-InsideMDAnderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/N_GTEx.tsv")

N_GTEx$Tissue[N_GTEx$Tissue=="Brain_Spinal_cord_cervical_c-1"] <- "Brain_Spinal_cord_cervical_c_1"
N_GTEx$Tissue[N_GTEx$Tissue=="Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV_transformed_lymphocytes"


method_consistency <- overlap_counts %>%
  group_by(Tissue, Method) %>%
  summarise(
    n_all4 = sum(n[overlap_group == "All 4 annotations"]),
    total_egenes = sum(n),
    prop_all4 = n_all4 / total_egenes,
    .groups = "drop"
  )

method_consistency <- method_consistency %>%
  mutate(
    discordance = 1 - prop_all4
  )

plotdat <- merge(x=method_consistency, y=N_GTEx,by="Tissue",all.x = T,all.y=F)


gg_1d_bottom <- ggplot(plotdat,
                       aes(x = N, y = discordance, color = Method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid")  +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),
    legend.position = "none",
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "Gene prediction discordance %\nacross annotation",
    color = "Method"
  ) +
  scale_y_continuous(labels = scales::percent_format())+
  scale_color_manual(values = c("featureCounts" = "#b9d47f",  # magenta
"kallisto"      = "#10b5a2", 
"Salmon"        = "#b5107b",  
"RSEM"          = "#6A3D9A"   
))



# figure 1a

gene_overlap <- long_genes %>%
  group_by(Tissue, Annotation, ensembl_gene_id) %>%
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
#   Annotation overlap_group      n  prop
#   <fct>      <fct>          <int> <dbl>
# 1 GENCODEv27 All 4 methods 321120 0.354
# 2 GENCODEv38 All 4 methods 331614 0.361
# 3 GENCODEv45 All 4 methods 335905 0.353
# 4 Ensembl    All 4 methods 341069 0.334


gg_1a_bottom <- ggplot(overlap_prop,
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
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f",  # magenta
"kallisto only"      = "#10b5a2",  
"Salmon only"        = "#b5107b",  
"RSEM only"          = "#6A3D9A",  
"All 4 methods" = "#80796b", 
"2 methods"    = "#56B4E9",   
"3 methods"    = "#0072B2"   
)) +
  labs(
    x = "Annotation", 
    y = expression("Proportion of genes with"~R^2 >= 0.01),
    fill = "Method"
  )

# figure 1c

N_GTEx <- read.delim("~/Library/CloudStorage/OneDrive-InsideMDAnderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/N_GTEx.tsv")

N_GTEx$Tissue[N_GTEx$Tissue=="Brain_Spinal_cord_cervical_c-1"] <- "Brain_Spinal_cord_cervical_c_1"
N_GTEx$Tissue[N_GTEx$Tissue=="Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV_transformed_lymphocytes"

method_consistency <- overlap_counts %>%
  group_by(Tissue, Annotation) %>%
  summarise(
    n_all4 = sum(n[overlap_group == "All 4 methods"]),
    total_egenes = sum(n),
    prop_all4 = n_all4 / total_egenes,
    .groups = "drop"
  )

method_consistency <- method_consistency %>%
  mutate(
    discordance = 1 - prop_all4
  )

plotdat <- merge(x=method_consistency, y=N_GTEx,by="Tissue",all.x = T,all.y=F)


gg_1c_bottom <- ggplot(plotdat,
                       aes(x = N, y = discordance, color = Annotation)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid")  +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),
    legend.position = "none",
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "Gene prediction discordance %\nacross method",
    color = "Annotation"
  ) +
  scale_y_continuous(labels = scales::percent_format())+
  scale_color_manual(values = c("Ensembl" = "#D55E00",   
"GENCODEv27"  = "#E69F00",   
"GENCODEv38"  = "#009E73",   
"GENCODEv45"  = "#CC79A7"   
))

# cowplot and print to PDF


fin_fig <- cowplot::plot_grid(
  gg_1a, gg_1b, gg_1c, gg_1d,
  gg_1a_bottom, gg_1b_bottom, gg_1c_bottom, gg_1d_bottom,
  nrow = 2,
  ncol = 4,
  rel_widths = c(1, 1, 0.6, 0.6),
  labels = c("A", "B", "C", "D", "", "", "", ""),
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1,
  label_size = 16,
  label_fontface = "bold"
)

ggsave(
  filename = "/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/Manuscript/Figures/r1/r1_sfig2_discordance_plot_all_annots_and_methods.pdf",
  plot = fin_fig,
  width = 11,    
  height = 9,     
  units = "in"
)

