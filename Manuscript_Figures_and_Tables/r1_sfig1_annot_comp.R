library(rtracklayer)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- 1. File paths ----------------------------------------------------------
gtf_dir  <- "/Users/sthead/OneDrive - Inside MD Anderson/gencode"
 
gtf_files <- list(
  "GENCODE v27"  = file.path(gtf_dir, "gencode.v27.annotation.gtf"),
  "GENCODE v38"  = file.path(gtf_dir, "gencode.v38.annotation.gtf"),
  "GENCODE v45"  = file.path(gtf_dir, "gencode.v45.annotation.gtf"),
  "Ensembl v115" = file.path(gtf_dir, "Ensembl_Homo_sapiens.GRCh38.115.chr.added.gtf")
)

read_gtf <- function(path, annot_label) {
  message("Reading: ", annot_label)
  df <- as.data.frame(import(path))
 
  ## Ensembl uses 'gene_biotype' / 'transcript_biotype'; GENCODE uses 'gene_type' / 'transcript_type'
  if (!"gene_type" %in% names(df) && "gene_biotype" %in% names(df)) {
    df$gene_type        <- df$gene_biotype
    df$transcript_type  <- df$transcript_biotype
  }
 
  ## Standardise NMD label  (Ensembl = "nonsense_mediated_decay")
  df$transcript_type <- gsub("nonsense_mediated_decay", "nonsense_mediated_decay", df$transcript_type)
 
  df$annotation <- annot_label
  df
}

raw_list <- mapply(read_gtf, gtf_files, names(gtf_files), SIMPLIFY = FALSE)

# Pull transcript IDs that have cds_start_NF or cds_end_NF directly from raw file
# these are partial transcripts

partial <- vector("list", 4)
names(partial) <- names(gtf_files)[1:4]

for(i in 1:4){
  gtf_path <- gtf_files[[i]]
  lines <- readLines(gtf_path)
  nf_lines <- lines[grepl("cds_start_NF|cds_end_NF", lines, fixed = FALSE)]
  nf_tx_ids <- unique(regmatches(nf_lines, 
                      regexpr('transcript_id "([^"]+)"', nf_lines)))
  nf_tx_ids <- gsub('transcript_id "|"', '', nf_tx_ids)
  partial[[i]] <- nf_tx_ids
}

## ---- 4. Summarise counts for Plot A -----------------------------------------
count_summary <- lapply(raw_list, function(df) {
  genes       <- df[df$type == "gene",       ]
  transcripts <- df[df$type == "transcript", ]
 
  lnc_genes <- sum(genes$gene_type %in% c("lncRNA", "lincRNA"),
                   na.rm = TRUE)
 
  lnc_tx    <- sum(transcripts$transcript_type %in%
                     c("lncRNA", "lincRNA"),
                   na.rm = TRUE)
 
  data.frame(
    annotation          = df$annotation[1],
    lncRNA_genes        = lnc_genes,
    lncRNA_transcripts  = lnc_tx,
    total_genes         = nrow(genes),
    total_transcripts   = nrow(transcripts)
  )
})
 
count_df <- bind_rows(count_summary)
 
## Long format for line plot
count_long <- count_df %>%
  pivot_longer(cols      = -annotation,
               names_to  = "metric",
               values_to = "count") %>%
  mutate(
    metric = recode(metric,
      lncRNA_genes       = "lncRNA genes",
      lncRNA_transcripts = "lncRNA transcripts",
      total_genes        = "Total genes",
      total_transcripts  = "Total transcripts"
    ),
    annotation = factor(annotation, levels = names(gtf_files))
  )
 
## ---- 5. Summarise transcript-type breakdown for Plot B ----------------------
  
tx_type_summary <- lapply(seq_along(raw_list[1:4]), function(i) {

  df <- raw_list[[i]]

  transcripts <- df[df$type == "transcript", ]

  NMD <- sum(

    transcripts$transcript_type == "nonsense_mediated_decay",

    na.rm = TRUE

  )

  PC <- unique(

    transcripts$transcript_id[

      transcripts$transcript_type == "protein_coding"

    ]

  )

  partial_length_PC <- sum(PC %in% partial[[i]])

  full_length_PC <- sum(!(PC %in% partial[[i]]))

  data.frame(

    annotation = names(partial)[i],

    full_length_PC = full_length_PC,

    NMD = NMD,

    partial_length_PC = partial_length_PC

  )

})

tx_type_df <- bind_rows(tx_type_summary)
 
tx_type_long <- tx_type_df %>%

  pivot_longer(

    cols = -annotation,

    names_to = "tx_type",

    values_to = "count"

  ) %>%

  mutate(

    tx_type = recode(

      tx_type,

      full_length_PC = "Full-length PC",

      NMD = "NMD",

      partial_length_PC = "Partial-length PC"

    ),

    annotation = factor(annotation, levels = names(gtf_files))

  )
 
## ---- 6. Compute ratios for Plot C ------------------------------------------
 
ratio_df <- lapply(raw_list, function(df) {
  genes       <- df[df$type == "gene",       ]
  transcripts <- df[df$type == "transcript", ]
 
  ## lncRNA biotypes
  lnc_biotypes <- c("lncRNA", "lincRNA")
 
  n_lnc_gene <- sum(genes$gene_type %in% lnc_biotypes, na.rm = TRUE)
  n_lnc_tx   <- sum(transcripts$transcript_type %in% lnc_biotypes |
                      (transcripts$gene_type %in% lnc_biotypes), na.rm = TRUE)
 
  n_pc_gene  <- sum(genes$gene_type == "protein_coding", na.rm = TRUE)
  n_pc_tx    <- sum(transcripts$transcript_type == "protein_coding", na.rm = TRUE)
 
  data.frame(
    annotation          = df$annotation[1],
    overall_ratio       = nrow(transcripts) / max(nrow(genes), 1),
    lncRNA_ratio        = n_lnc_tx / max(n_lnc_gene, 1),
    PC_ratio            = n_pc_tx  / max(n_pc_gene,  1)
  )
}) %>% bind_rows()
 
ratio_long <- ratio_df %>%
  pivot_longer(cols = -annotation, names_to = "ratio_type", values_to = "ratio") %>%
  mutate(
    ratio_type = recode(ratio_type,
      overall_ratio = "Transcripts per gene",
      lncRNA_ratio  = "lncRNA transcripts per lncRNA gene",
      PC_ratio      = "PC transcripts per PC gene"
    ),
    annotation = factor(annotation, levels = names(gtf_files))
  )
 
## ---- 7. Plot A: Line plot ---------------------------------------------------
 
## Colour palette (colour-blind-friendly)
line_colours <- c(
  "lncRNA genes"        = "#3FB0E0",
  "lncRNA transcripts"  = "#7F3FE0",
  "Total genes"         = "#48BD7F",
  "Total transcripts"   = "#D68927"
)
line_shapes <- c(
  "lncRNA genes"        = 16,
  "lncRNA transcripts"  = 17,
  "Total genes"         = 15,
  "Total transcripts"   = 18
)
 
pA <- ggplot(count_long,
             aes(x = annotation, y = count, colour = metric,
                 group = metric, shape = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = line_colours) +
  scale_shape_manual(values  = line_shapes) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title   = "Growth trends across annotations",
    x       = NULL,
    y       = "Count",
    colour  = NULL,
    shape   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 13),
  )
 
## ---- 8. Plot B: Stacked bar — transcript types per annotation ---------------
 
bar_colours_B <- c(
  "Full-length PC"    = "#509E0B",
  "NMD"               = "#E0971F",
  "Partial-length PC" = "#E6DB6A"
)
 
pB <- ggplot(tx_type_long,
             aes(x = annotation, y = count, fill = tx_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = bar_colours_B) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Protein-coding and NMD transcripts",
    x     = NULL,
    y     = "Number of transcripts",
    fill  = "Transcript type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )
 
## ---- 9. Plot C: Bar plot — transcript-to-gene ratios -----------------------
 
 
pC <- ggplot(ratio_long,
             aes(x = annotation, y = ratio, fill = ratio_type)) +
  geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
  labs(
    title = "Transcript-to-gene ratios by annotation",
    x     = NULL,
    y     = "Transcripts per gene",
    fill  = "Ratio type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )
 
## ---- 10. Combine and save ---------------------------------------------------


fin_fig <- cowplot::plot_grid(
  pA, pB,pC,
  nrow = 3,
  ncol = 1,
  rel_widths = c(1, 1, 1),
  labels = c("A", "B", "C"),
  label_x = 0.01,
  label_y = 1,
  hjust = 0,
  vjust = 1,
  label_size = 16,
  label_fontface = "bold"
)

ggsave(
  filename = "/Users/sthead/OneDrive - Inside MD Anderson/Bhattacharya,Arjun's files - GTEx GENCODE Comp/Manuscript/Figures/r1/r1_sfig1_annot_comp.pdf",
  plot = fin_fig,
  width = 8,     # adjust to taste
  height = 9,     # adjust to taste
  units = "in"
)
 
 

