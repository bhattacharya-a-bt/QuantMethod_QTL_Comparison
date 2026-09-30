library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggtranscript)
library(pheatmap)
library(GenomicRanges)
library(gridExtra)
library(dplyr)
library(viridis)
library(reshape2)
library(ggdendro)
library(cowplot)
library(ggnewscale)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(ggpubr)
library(forcats)

gtf_gencode_v45 <- import("gencode.v45.annotation.gtf") # download file from https://www.gencodegenes.org/human/releases.html
gtf_gencode_v38 <- import("gencode.v38.annotation.gtf") # download file from https://www.gencodegenes.org/human/releases.html
gtf_gencode_v27 <- import("gencode.v27.annotation.gtf") # download file from https://www.gencodegenes.org/human/releases.html

extract_gene_tx_features <- function(gtf_gr, reference_name) {
  
  # 1) selcect gene / transcript / exon
  genes <- gtf_gr[gtf_gr$type == "gene"]
  tx    <- gtf_gr[gtf_gr$type == "transcript"]
  exons <- gtf_gr[gtf_gr$type == "exon"]
  
  # 2) gene-level information
  gene_df <- as.data.frame(genes) %>%
    transmute(
      gene_id   = gene_id,
      gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
      biotype   = gene_type,
      chr       = as.character(seqnames),
      start     = start,
      end       = end,
      gene_length = end - start + 1
    ) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  # 3) transcript
  tx_df <- as.data.frame(tx) %>%
    select(gene_id, transcript_id, width) %>%
    distinct()
  
  tx_count <- tx_df %>%
    count(gene_id, name = "transcript_count")
  
  # 4) exon per gene
  exon_df <- as.data.frame(exons) %>%
    select(gene_id, exon_id) %>%
    distinct()
  
  exon_count <- exon_df %>%
    count(gene_id, name = "exon_count")
  
  
  gene_annot <- gene_df %>%
    left_join(tx_count,  by = "gene_id") %>%
    left_join(exon_count, by = "gene_id") %>%
    mutate(
      transcript_count = tidyr::replace_na(transcript_count, 0L),
      exon_count       = tidyr::replace_na(exon_count, 0L),
      
      is_noncoding     = ifelse(str_detect(biotype, "protein_coding"), 0L, 1L),
      reference        = reference_name
    )
  
  return(gene_annot)
}

gene_v27 <- extract_gene_tx_features(gtf_gencode_v27, "v27")
gene_v38 <- extract_gene_tx_features(gtf_gencode_v38, "v38")
gene_v45 <- extract_gene_tx_features(gtf_gencode_v45, "v45")

genes_all <- bind_rows(gene_v27, gene_v38, gene_v45)
rm(gtf_gencode_v27)
rm(gtf_gencode_v38)
rm(gtf_gencode_v45)


###### N exons and gene length

# across methods

egene <- readRDS("figure_r1_aggregated_eGene_lists.RDS")
total <- egene %>%
  unnest(ensembl_gene_id) %>%
  filter(Annotation != "Ensembl")

egene_class <- total %>%
  select(Method, Annotation, ensembl_gene_id) %>%
  distinct() %>%
  mutate(
    Annotation = case_when(
      Annotation == "GENCODE_v27" ~ "v27",
      Annotation == "GENCODE_v38" ~ "v38",
      Annotation == "GENCODE_v45" ~ "v45",
      TRUE ~ Annotation
    )
    
  ) %>%
  group_by(Method, ensembl_gene_id) %>%
  summarise(
    n_Annotations = n_distinct(Annotation),
    Annotations   = paste(sort(unique(Annotation)), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    n_Annotations == 4 ~ "shared_4_Annotations",
    n_Annotations == 3 ~ "shared_3_Annotations",
    n_Annotations == 2 ~ "shared_2_Annotations",
    n_Annotations == 1 ~ "unique_1_Annotation"
  ))


egene_vs <- egene_class %>%
  mutate(
    Annotation_set = Annotations %>%
      str_split(",") %>% lapply(sort) %>% sapply(paste, collapse = ",")
  )

## ---- 1) Clean Annotation_set formatting  ----
egene_vs <- egene_vs %>%
  mutate(
    Annotation_set = str_replace_all(Annotation_set, "\\s+", ""),
  )

Annotation_levels <- c("v27", "v38", "v45", "v27,v38", "v27,v45", "v38,v45", "v27,v38,v45")

method_levels <- c("kallisto","salmon","featureCounts","RSEM")


## ---- 1) build egene_annot for reference ----
make_egene_annot_by_reference <- function(egene_vs, genes_all_ref,
                                          Annotation_levels, method_levels) {
  
  genes_ref <- genes_all_ref %>%
    mutate(ensembl_gene_id = str_remove(gene_id, "\\..*$")) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  egene_annot <- egene_vs %>%
    mutate(
      Annotation_set = factor(Annotation_set, levels = Annotation_levels),
      Method = factor(Method, levels = method_levels)
    ) %>%
    left_join(genes_ref, by = "ensembl_gene_id")
  
  #message("Missing gene_length after join: ", sum(is.na(egene_annot$gene_length)))
  egene_annot
}

plot_structure_metric <- function(egene_annot, metric,
                                  ref_label,
                                  log10_y = FALSE,
                                  y_label = NULL,
                                  title = NULL) {
  
  stopifnot(metric %in% c("gene_length", "transcript_count", "exon_count"))
  
  # count_labels <- egene_annot |>
  # dplyr::group_by(Method, Annotation_set) |>
  # dplyr::summarise(
  #   n = sum(!is.na(.data[[metric]])),
  #   y_min = 8,
  #   y_max = max(.data[[metric]], na.rm = TRUE),
  #   .groups = "drop"
  # ) |>
  # dplyr::mutate(
  #   label_y = if (log10_y) {
  #     y_min / 1.5
  #   } else {
  #     y_min - 0.08 * (y_max - y_min)
  #   }
  # )

  p <- ggplot(
    egene_annot,
    aes(x = Annotation_set, y = .data[[metric]], fill = Method)
  ) +
    geom_boxplot(
      outlier.size = 0.4,
      alpha = 0.7,
      width = 0.7
    ) +
#     geom_text(
#   data = count_labels,
#   aes(x = Annotation_set, y = label_y, label = paste0(" n=", n)),
#   inherit.aes = FALSE,
#   size = 2, angle=90,vjust=0.5
# ) +
    facet_wrap(~ Method, nrow = 1) +
    scale_fill_manual(
      values = c(
        "salmon" = "#1b9e77",
        "RSEM"  = "#d95f02",
        "kallisto" = "#0099cc"      ),
      guide = "none"   # legend unnecessary because of facets
    ) +
    labs(
      x = "GENCODE annotation",
      y = ifelse(is.null(y_label), metric, y_label),
      title = paste0(
        title,
        "Gene annotation information: ",
        ref_label
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  if (log10_y) {
    p <- p + scale_y_log10(labels = scales::comma)
  }
  
  p
}



`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- 3) Split references ----

genes_all_v45 <- subset(genes_all, reference == "v45")
genes_all_v38 <- subset(genes_all, reference == "v38")
genes_all_v27 <- subset(genes_all, reference == "v27")

## ---- 4) Build egene_annot for each reference ----

egene_annot_v45 <- make_egene_annot_by_reference(egene_vs, genes_all_v45, Annotation_levels, method_levels)
egene_annot_v38 <- make_egene_annot_by_reference(egene_vs, genes_all_v38, Annotation_levels, method_levels)
egene_annot_v27 <- make_egene_annot_by_reference(egene_vs, genes_all_v27, Annotation_levels, method_levels)

## ---- 5) (A) Gene length plots (log10 scale recommended) ----

p_len_v27 <- plot_structure_metric(egene_annot_v27, "gene_length", "GENCODE v27",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")
p_len_v38 <- plot_structure_metric(egene_annot_v38, "gene_length", "GENCODE v38",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")
p_len_v45 <- plot_structure_metric(egene_annot_v45, "gene_length", "GENCODE v45",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")

fig_len <- ggarrange(p_len_v27, p_len_v38, p_len_v45,
                     ncol = 1, nrow = 3, align = "v")
ggsave("r1_sfig10_discord_by_genelength_ax_annot.pdf", fig_len, width = 10, height = 14)


## ---- 6) (B) Transcript count plots ----

p_tx_v27 <- plot_structure_metric(egene_annot_v27, "transcript_count", "GENCODE v27",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")
p_tx_v38 <- plot_structure_metric(egene_annot_v38, "transcript_count", "GENCODE v38",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")
p_tx_v45 <- plot_structure_metric(egene_annot_v45, "transcript_count", "GENCODE v45",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")

fig_tx <- ggarrange(p_tx_v27, p_tx_v38, p_tx_v45,
                    ncol = 1, nrow = 3, align = "v",common.legend = TRUE)

ggsave("r1_sfig9_discord_by_txct_ax_annot.pdf", fig_tx, width = 10, height = 14)

## ---- 7) (C) Exon count plots ----

p_ex_v27 <- plot_structure_metric(egene_annot_v27, "exon_count", "GENCODE v27",
                                  log10_y = TRUE, y_label = "Exon count (log10)")
p_ex_v38 <- plot_structure_metric(egene_annot_v38, "exon_count", "GENCODE v38",
                                  log10_y = TRUE, y_label = "Exon count (log10)")
p_ex_v45 <- plot_structure_metric(egene_annot_v45, "exon_count", "GENCODE v45",
                                  log10_y = TRUE, y_label = "Exon count (log10)")

fig_ex <- ggarrange(p_ex_v27, p_ex_v38, p_ex_v45,
                    ncol = 1, nrow = 3, align = "v")

ggsave("r1_sfig8_discord_by_exonct_ax_annot.pdf", fig_ex, width = 10, height = 14)





# across annotations

egene <- readRDS("figure_r1_aggregated_eGene_lists.RDS")
total <- egene %>%
  unnest(ensembl_gene_id) %>%
  filter(Annotation != "Ensembl")

egene_class <- total %>%
  select(Method, Annotation, ensembl_gene_id) %>%
  distinct() %>%
  mutate(
    Annotation = case_when(
      Annotation == "GENCODE_v27" ~ "v27",
      Annotation == "GENCODE_v38" ~ "v38",
      Annotation == "GENCODE_v45" ~ "v45",
      TRUE ~ Annotation
    )
    
  ) %>%
  group_by(Annotation, ensembl_gene_id) %>%
  summarise(
    n_Annotations = n_distinct(Method),
    Annotations   = paste(sort(unique(Method)), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    n_Annotations == 4 ~ "shared_4_Annotations",
    n_Annotations == 3 ~ "shared_3_Annotations",
    n_Annotations == 2 ~ "shared_2_Annotations",
    n_Annotations == 1 ~ "unique_1_Annotation"
  ))


egene_vs <- egene_class %>%
  mutate(
    Annotation_set = Annotations %>%
      str_split(",") %>% lapply(sort) %>% sapply(paste, collapse = ",")
  )

## ---- 1) Clean Annotation_set formatting  ----

egene_vs <- egene_vs %>%

  mutate(

    comma_count = str_count(Annotation_set, ",")

  ) %>%

  mutate(

    Annotation_set = fct_relevel(

      Annotation_set,

      unique(Annotation_set[order(comma_count)])

    )

  ) %>%

  select(-comma_count)  

Annotation_levels <- levels(egene_vs$Annotation_set)

method_levels <- c("v27","v38","v45")


## ---- 1) build egene_annot for reference ----
make_egene_annot_by_reference <- function(egene_vs, genes_all_ref,
                                          Annotation_levels, method_levels) {
  
  genes_ref <- genes_all_ref %>%
    mutate(ensembl_gene_id = str_remove(gene_id, "\\..*$")) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  egene_annot <- egene_vs %>%
    mutate(
      Annotation_set = factor(Annotation_set, levels = Annotation_levels),
      Method = factor(Annotation, levels = method_levels)
    ) %>%
    left_join(genes_ref, by = "ensembl_gene_id")
  
  #message("Missing gene_length after join: ", sum(is.na(egene_annot$gene_length)))
  egene_annot
}

plot_structure_metric <- function(egene_annot, metric,
                                  ref_label,
                                  log10_y = FALSE,
                                  y_label = NULL,
                                  title = NULL) {
  
  stopifnot(metric %in% c("gene_length", "transcript_count", "exon_count"))
  
  p <- ggplot(
    egene_annot,
    aes(x = Annotation_set, y = .data[[metric]], fill = Method)
  ) +
    geom_boxplot(
      outlier.size = 0.4,
      alpha = 0.7,
      width = 0.7
    ) +
    facet_wrap(~ Method, nrow = 1) +
    scale_fill_manual(
      values = c(
        "v27" = "#1b9e77",
        "v38"  = "#d95f02",
        "v45" = "#0099cc"      ),
      guide = "none"   # legend unnecessary because of facets
    ) +
    labs(
      x = "Method",
      y = ifelse(is.null(y_label), metric, y_label),
      title = paste0(
        title,
        "Gene annotation information: ",
        ref_label
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  if (log10_y) {
    p <- p + scale_y_log10(labels = scales::comma)
  }
  
  p
}


`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- 3) Split references ----

genes_all_v45 <- subset(genes_all, reference == "v45")
genes_all_v38 <- subset(genes_all, reference == "v38")
genes_all_v27 <- subset(genes_all, reference == "v27")

## ---- 4) Build egene_annot for each reference ----

egene_annot_v45 <- make_egene_annot_by_reference(egene_vs, genes_all_v45, Annotation_levels, method_levels)
egene_annot_v38 <- make_egene_annot_by_reference(egene_vs, genes_all_v38, Annotation_levels, method_levels)
egene_annot_v27 <- make_egene_annot_by_reference(egene_vs, genes_all_v27, Annotation_levels, method_levels)

## ---- 5) (A) Gene length plots (log10 scale recommended) ----

p_len_v27 <- plot_structure_metric(egene_annot_v27, "gene_length", "GENCODE v27",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")
p_len_v38 <- plot_structure_metric(egene_annot_v38, "gene_length", "GENCODE v38",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")
p_len_v45 <- plot_structure_metric(egene_annot_v45, "gene_length", "GENCODE v45",
                                   log10_y = TRUE, y_label = "Gene length (bp, log10)")

fig_len <- ggarrange(p_len_v27, p_len_v38, p_len_v45,
                     ncol = 1, nrow = 3, align = "v")
ggsave("r1_sfig13_discord_by_genelength_ax_methods.pdf", fig_len, width = 10, height = 14)


## ---- 6) (B) Transcript count plots ----

p_tx_v27 <- plot_structure_metric(egene_annot_v27, "transcript_count", "GENCODE v27",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")
p_tx_v38 <- plot_structure_metric(egene_annot_v38, "transcript_count", "GENCODE v38",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")
p_tx_v45 <- plot_structure_metric(egene_annot_v45, "transcript_count", "GENCODE v45",
                                  log10_y = TRUE, y_label = "Transcript count (log10)")

fig_tx <- ggarrange(p_tx_v27, p_tx_v38, p_tx_v45,
                    ncol = 1, nrow = 3, align = "v",common.legend = TRUE)

ggsave("r1_sfig12_discord_by_txct_ax_methods.pdf", fig_tx, width = 10, height = 14)

## ---- 7) (C) Exon count plots ----

p_ex_v27 <- plot_structure_metric(egene_annot_v27, "exon_count", "GENCODE v27",
                                  log10_y = TRUE, y_label = "Exon count (log10)")
p_ex_v38 <- plot_structure_metric(egene_annot_v38, "exon_count", "GENCODE v38",
                                  log10_y = TRUE, y_label = "Exon count (log10)")
p_ex_v45 <- plot_structure_metric(egene_annot_v45, "exon_count", "GENCODE v45",
                                  log10_y = TRUE, y_label = "Exon count (log10)")

fig_ex <- ggarrange(p_ex_v27, p_ex_v38, p_ex_v45,
                    ncol = 1, nrow = 3, align = "v")

ggsave("r1_sfig11_discord_by_exonct_ax_methods.pdf", fig_ex, width = 10, height = 14)


