library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(data.table)
library(cowplot)

# load for for top row 

dat_gtex <- readRDS("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_aggregated_eGene_lists.RDS")
#dat_gtex <- dat_gtex[dat_gtex$Tissue=="Breast_Mammary_Tissue",]
dat_gtex <- dat_gtex[dat_gtex$Tissue %in% c("Breast_Mammary_Tissue","Uterus","Vagina"),]
dat_gtex$Source <- "GTEx Breast"
dat_gtex$Source[dat_gtex$Tissue=="Uterus"] <- "GTEx Uterus"
dat_gtex$Source[dat_gtex$Tissue=="Vagina"] <- "GTEx Vagina"
dat_gtex <- dat_gtex[,-which(names(dat_gtex)=="Tissue")]

dat <- readRDS("/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/3_results/r1_wz_aggregated_eGene_lists.RDS")
dat$Source <- "Ping et al. Breast"

dat_sub <- rbind(dat,dat_gtex)
dat_sub <- dat_sub[grep("GENCODE",dat_sub$Annotation),]
dat_sub <- dat_sub[dat_sub$Method %in% c("featureCounts","salmon"),]

dat_sub$Annotation <- factor(dat_sub$Annotation,levels=c("GENCODE_v27","GENCODE_v38","GENCODE_v45"),
                             labels=c("v27","v38","v45"))
dat_sub$Method <- factor(dat_sub$Method,levels=c("salmon","featureCounts"),
                         labels=c("Salmon","featureCounts"))
dat_sub$Source <- factor(dat_sub$Source,levels=c("Ping et al. Breast","GTEx Breast","GTEx Uterus","GTEx Vagina"),
  labels=c("Ping et al. Breast (N=150)","GTEx Breast (N=459)","GTEx Uterus (N=142)","GTEx Vagina (N=156)"))

long_genes <- dat_sub %>%
  select(Source, Annotation, Method, ensembl_gene_id) %>%
  unnest(ensembl_gene_id)

# figure 1b

gene_overlap<- long_genes %>%
  group_by(Source, Method, ensembl_gene_id) %>%
  summarise(
    annotations = list(sort(unique(Annotation))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(annotations, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 3) return("All three")
      paste(x, collapse = " & ")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Source, Method, overlap_group)

overlap_counts <- overlap_counts %>%
  mutate(
    overlap_group = factor(
      overlap_group,
      levels = c(
        # one annotation only
        "v27 only",
        "v38 only",
        "v45 only",
        
        # two annotations
        "v27 & v38",
        "v27 & v45",
        "v38 & v45",
        
        # four annotations
        "All three"
      )
    )
  )

overlap_counts$overlap_group <- factor(overlap_counts$overlap_group,levels=c(
    "v27 only","v38 only", "v45 only","v27 & v38",
        "v27 & v45",
        "v38 & v45","All three"
    ))


overlap_prop <- overlap_counts %>%
  group_by(Source,Method, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Source,Method) %>%
  mutate(prop = n / sum(n))

gg_1b <- ggplot(overlap_prop,
                aes(x = Method, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") +facet_grid(.~Source)+
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
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("All three" = "#80796b", "v27 only" = "#E69F00", "v38 only" = "#009E73", "v45 only"="#CC79A7",
                               "v27 & v38"="#b24745", "v27 & v45"="#79af97", "v38 & v45"="#6a6599")) +
  labs(
    x = "Method", 
    y = "Proportion of eGenes",
    fill = "GENCODE"
  )


# figure 1a

gene_overlap <- long_genes %>%
  group_by(Source, Annotation, ensembl_gene_id) %>%
  summarise(
    methods = list(sort(unique(Method))),
    .groups = "drop"
  ) %>%
  mutate(
    # extract overlap_group as a single string
    overlap_group = sapply(methods, function(x) {
      if (length(x) == 1) return(paste0(x, " only"))
      if (length(x) == 2) return("Both")
      paste(x, collapse = " & ")
    })
  )

overlap_counts <- gene_overlap %>%
  count(Source, Annotation, overlap_group)

overlap_counts <- overlap_counts %>%
  mutate(
    overlap_group = factor(
      overlap_group,
      levels = c(
        # one annotation only
        "featureCounts only",
        "Salmon only",
        
        # four annotations
        "Both"
      )
    )
  )


overlap_prop <- overlap_counts %>%
  group_by(Source,Annotation, overlap_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Source,Annotation) %>%
  mutate(prop = n / sum(n))

gg_1a <- ggplot(overlap_prop,
                aes(x = Annotation, y = prop, fill = overlap_group)) +
  geom_bar(stat = "identity") + facet_grid(.~Source)+
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
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("featureCounts only" = "#b9d47f", "Salmon only" = "#b5107b", "Both" = "#80796b")) +
  labs(
    x = "GENCODE", 
    y = "Proportion of eGenes",
    fill = "Method"
  )


# cowplot and print to PDF


fin_fig <- cowplot::plot_grid(
  gg_1b, gg_1a, 
  nrow = 2,
  ncol = 1,
  rel_widths=c(1,1)
)

fin_fig <- cowplot::plot_grid(
  gg_1b, 
  gg_1a,
  nrow = 2,
  ncol = 1,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_x = 0.0,
  label_y = .995,
  hjust = 0,
  vjust = 1,
  label_size = 16,
  label_fontface = "bold")


ggsave(
  filename = "/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/Manuscript/Figures/r1/r1_sfig_wzhang.pdf",
  plot = fin_fig,
  width = 10,     # adjust to taste
  height = 10,     # adjust to taste
  units = "in"
)


