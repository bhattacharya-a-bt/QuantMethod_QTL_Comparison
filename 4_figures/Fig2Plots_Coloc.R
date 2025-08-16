setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/2_qtlanalysis/')

### Load libraries
library(dplyr)
library(ggplot2)
library(grid)
library(tidyr)
library(paletteer)
library(scales)
library(gridExtra)
library(hexbin)  # For hexagonal binning
library(cowplot)
library(data.table)

coloc_res = fread('colocalization_results_combined.tsv')


# Parse chromosome and position from GWAS_CHRPOS
coloc_res <- coloc_res %>%
  separate(GWAS_CHRPOS, into = c("chr_name", "position"), sep = ":", convert = TRUE) %>%
  mutate(
    chr_name = as.numeric(chr_name),
    position = as.numeric(position)
  ) %>%
  filter(!is.na(chr_name) & !is.na(position))  # Remove rows with parsing issues

# Create cancer name mapping
cancer_names <- c(
  "Colorectal" = "Colorectal",
  "ProstateOverall" = "Prostate", 
  "BreastOverall" = "Breast",
  "LungOverall" = "Lung",
  "OvarianOverall" = "Ovarian",
  "Melanoma" = "Melanoma",
  "Glioma" = "Glioma",
  "Renal" = "Renal",
  "HeadNeck" = "HeadNeck",
  "Endometrial" = "Endometrial"
)

# Define GWAS loci (SNPs within 1MB of each other) grouped by cancer
coloc_res_with_loci <- coloc_res %>%
  arrange(Cancer, chr_name, position) %>%
  group_by(Cancer) %>%
  mutate(
    # Create locus groups: new locus if >1MB from previous SNP or different chromosome
    position_diff = c(1e6 + 1, diff(position)),  # First SNP gets large diff to start new locus
    chr_diff = c(1, diff(chr_name)),
    new_locus = (position_diff > 1e6) | (chr_diff != 0),
    locus_id = cumsum(new_locus),
    gwas_locus = paste(Cancer, chr_name, locus_id, sep = "_")
  ) %>%
  ungroup()

# Create method overlap analysis for genes
method_overlap_coloc <- coloc_res_with_loci %>%
  group_by(Cancer, Version, Gene) %>%
  summarise(
    methods = paste(sort(unique(Method)), collapse = ","),
    n_methods = n_distinct(Method),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_group = case_when(
      methods == "salmon" ~ "Salmon only",
      methods == "star" ~ "STAR only", 
      methods == "salmon,star" ~ "Both",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Cancer, Version, overlap_group) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(Cancer, Version) %>%
  mutate(
    total_genes = sum(count),
    proportion = count / total_genes
  ) %>%
  ungroup()

# Create version overlap analysis for genes
version_overlap_coloc <- coloc_res_with_loci %>%
  group_by(Cancer, Method, Gene) %>%
  summarise(
    versions = paste(sort(unique(Version)), collapse = ","),
    n_versions = n_distinct(Version),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_group = case_when(
      versions == "v27" ~ "v27 only",
      versions == "v38" ~ "v38 only",
      versions == "v45" ~ "v45 only",
      versions == "v27,v38" ~ "v27 & v38",
      versions == "v27,v45" ~ "v27 & v45", 
      versions == "v38,v45" ~ "v38 & v45",
      versions == "v27,v38,v45" ~ "All three",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Cancer, Method, overlap_group) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(Cancer, Method) %>%
  mutate(
    total_genes = sum(count),
    proportion = count / total_genes
  ) %>%
  ungroup()

# Prepare data for plotting - genes
method_plot_data_coloc <- method_overlap_coloc %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", gsub("v", "", Version)),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

version_plot_data_coloc <- version_overlap_coloc %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = toupper(Method),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

# Create method overlap analysis for GWAS loci
method_overlap_loci_coloc <- coloc_res_with_loci %>%
  group_by(Cancer, Version, gwas_locus) %>%
  summarise(
    methods = paste(sort(unique(Method)), collapse = ","),
    n_methods = n_distinct(Method),
    chr_name = first(chr_name),
    position = first(position),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_group = case_when(
      methods == "salmon" ~ "Salmon only",
      methods == "star" ~ "STAR only", 
      methods == "salmon,star" ~ "Both",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Cancer, Version, overlap_group) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(Cancer, Version) %>%
  mutate(
    total_loci = sum(count),
    proportion = count / total_loci
  ) %>%
  ungroup()

# Create version overlap analysis for GWAS loci
version_overlap_loci_coloc <- coloc_res_with_loci %>%
  group_by(Cancer, Method, gwas_locus) %>%
  summarise(
    versions = paste(sort(unique(Version)), collapse = ","),
    n_versions = n_distinct(Version),
    chr_name = first(chr_name),
    position = first(position),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_group = case_when(
      versions == "v27" ~ "v27 only",
      versions == "v38" ~ "v38 only",
      versions == "v45" ~ "v45 only",
      versions == "v27,v38" ~ "v27 & v38",
      versions == "v27,v45" ~ "v27 & v45", 
      versions == "v38,v45" ~ "v38 & v45",
      versions == "v27,v38,v45" ~ "All three",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Cancer, Method, overlap_group) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(Cancer, Method) %>%
  mutate(
    total_loci = sum(count),
    proportion = count / total_loci
  ) %>%
  ungroup()

# Prepare data for plotting - GWAS loci
method_plot_data_loci_coloc <- method_overlap_loci_coloc %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", gsub("v", "", Version)),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

version_plot_data_loci_coloc <- version_overlap_loci_coloc %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = toupper(Method),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

# Summary plots - aggregate across cancers - genes
method_summary_data_coloc <- method_plot_data_coloc %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

version_summary_data_coloc <- version_plot_data_coloc %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

# Summary plots - aggregate across cancers - GWAS loci
method_summary_data_loci_coloc <- method_plot_data_loci_coloc %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

version_summary_data_loci_coloc <- version_plot_data_loci_coloc %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

# GENE PLOTS
# Plot 1: Method comparison by cancer - stacked counts (genes)
p1_cancer_genes <- ggplot(method_plot_data_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Version, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Number of Colocalized Genes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison by cancer - proportions (genes)
p2_cancer_genes <- ggplot(method_plot_data_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Version, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Proportion of Colocalized Genes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison by cancer - stacked counts (genes)
p3_cancer_genes <- ggplot(version_plot_data_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Number of Colocalized Genes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison by cancer - proportions (genes)
p4_cancer_genes <- ggplot(version_plot_data_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Method, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Proportion of Colocalized Genes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots (main figures) - genes
p1_summary_genes <- ggplot(method_summary_data_coloc, aes(x = Version, y = total_count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
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
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0,0,0,0)
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "GENCODE", 
    y = "Proportion of Colocalized Genes",
    fill = "Method"
  )

p2_summary_genes <- ggplot(version_summary_data_coloc, aes(x = Method, y = total_count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
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
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0,0,0,0)
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Method",
    y = "Proportion of Colocalized Genes", 
    fill = "GENCODE"
  )

# GWAS LOCI PLOTS
# Plot 1: Method comparison by cancer - stacked counts (loci)
p1_cancer_loci <- ggplot(method_plot_data_loci_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Version, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Number of GWAS Loci",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison by cancer - proportions (loci)
p2_cancer_loci <- ggplot(method_plot_data_loci_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Version, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Proportion of GWAS Loci\ntagged by colocalization",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison by cancer - stacked counts (loci)
p3_cancer_loci <- ggplot(version_plot_data_loci_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Number of GWAS Loci",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison by cancer - proportions (loci)
p4_cancer_loci <- ggplot(version_plot_data_loci_coloc, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Method, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Proportion of GWAS Loci\ntagged by colocalization",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots (main figures) - loci
p1_summary_loci <- ggplot(method_summary_data_loci_coloc, aes(x = Version, y = total_count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
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
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0,0,0,0)
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "GENCODE", 
    y = "Proportion of GWAS Loci\ntagged by colocalization",
    fill = "Method"
  )

p2_summary_loci <- ggplot(version_summary_data_loci_coloc, aes(x = Method, y = total_count, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "fill") +
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
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0,0,0,0)
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Method",
    y = "Proportion of GWAS Loci\ntagged by colocalization", 
    fill = "GENCODE"
  )

# Create combined main figures
p_genes_main <- plot_grid(p1_summary_genes, p2_summary_genes, ncol = 2)
p_loci_main <- plot_grid(p1_summary_loci, p2_summary_loci, ncol = 2)

# Print summary statistics
cat("\n=== COLOCALIZATION SUMMARY ===\n")
total_genes_by_cancer <- coloc_res_with_loci %>%
  group_by(Cancer) %>%
  summarise(unique_genes = n_distinct(Gene), .groups = "drop") %>%
  mutate(Cancer_clean = cancer_names[Cancer])

total_loci_by_cancer <- coloc_res_with_loci %>%
  group_by(Cancer) %>%
  summarise(unique_loci = n_distinct(gwas_locus), .groups = "drop") %>%
  mutate(Cancer_clean = cancer_names[Cancer])

print("Genes by cancer:")
print(total_genes_by_cancer)

print("GWAS loci by cancer:")
print(total_loci_by_cancer)

cat(paste("\nTotal unique colocalized genes across all cancers:", 
          n_distinct(coloc_res_with_loci$Gene), "\n"))
cat(paste("Total unique GWAS loci across all cancers:", 
          n_distinct(coloc_res_with_loci$gwas_locus), "\n"))

# Save plots to specified directory
output_dir <- "/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/"

# Save gene plots
ggsave(plot = p_genes_main,
       filename = paste0(output_dir, "Coloc_Genes_Method_Version_Comparison_Main.pdf"),
       height = 4,
       width = 7)

ggsave(plot = p1_cancer_genes,
       filename = paste0(output_dir, "Coloc_Genes_Method_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p2_cancer_genes,
       filename = paste0(output_dir, "Coloc_Genes_Method_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p3_cancer_genes,
       filename = paste0(output_dir, "Coloc_Genes_Version_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p4_cancer_genes,
       filename = paste0(output_dir, "Coloc_Genes_Version_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

# Save GWAS loci plots
ggsave(plot = p_loci_main,
       filename = paste0(output_dir, "Coloc_GWAS_Loci_Method_Version_Comparison_Main.pdf"),
       height = 4,
       width = 7)

ggsave(plot = p1_cancer_loci,
       filename = paste0(output_dir, "Coloc_GWAS_Loci_Method_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p2_cancer_loci,
       filename = paste0(output_dir, "Coloc_GWAS_Loci_Method_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p3_cancer_loci,
       filename = paste0(output_dir, "Coloc_GWAS_Loci_Version_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p4_cancer_loci,
       filename = paste0(output_dir, "Coloc_GWAS_Loci_Version_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

# Export tables to CSV files
output_dir_tables <- "/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/"

# Table 1: Number of colocalized genes per cancer - Method overlap
coloc_genes_method_table <- method_plot_data_coloc %>%
  dplyr::select(Cancer_clean, Version, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean,
    `GENCODE Version` = Version
  ) %>%
  arrange(`Cancer Type`, `GENCODE Version`)

# Table 2: Number of colocalized genes per cancer - Version overlap
coloc_genes_version_table <- version_plot_data_coloc %>%
  dplyr::select(Cancer_clean, Method, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean
  ) %>%
  arrange(`Cancer Type`, Method)

# Table 3: Number of GWAS loci per cancer - Method overlap
coloc_loci_method_table <- method_plot_data_loci_coloc %>%
  dplyr::select(Cancer_clean, Version, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean,
    `GENCODE Version` = Version
  ) %>%
  arrange(`Cancer Type`, `GENCODE Version`)

# Table 4: Number of GWAS loci per cancer - Version overlap
coloc_loci_version_table <- version_plot_data_loci_coloc %>%
  dplyr::select(Cancer_clean, Method, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean
  ) %>%
  arrange(`Cancer Type`, Method)

# Write tables to CSV files
write.csv(coloc_genes_method_table, 
          paste0(output_dir_tables, "Table_Coloc_Genes_Method_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(coloc_genes_version_table,
          paste0(output_dir_tables, "Table_Coloc_Genes_Version_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(coloc_loci_method_table,
          paste0(output_dir_tables, "Table_Coloc_GWAS_Loci_Method_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(coloc_loci_version_table,
          paste0(output_dir_tables, "Table_Coloc_GWAS_Loci_Version_Overlap_by_Cancer.csv"),
          row.names = FALSE)

# Print summary of what was saved
cat("\n=== PLOTS AND TABLES EXPORTED ===\n")
cat("GENE PLOTS:\n")
cat("1. Coloc_Genes_Method_Version_Comparison_Main.pdf - Combined summary plot\n")
cat("2. Coloc_Genes_Method_Comparison_by_Cancer_Counts.pdf - Method comparison by cancer (counts)\n")
cat("3. Coloc_Genes_Method_Comparison_by_Cancer_Proportions.pdf - Method comparison by cancer (proportions)\n")
cat("4. Coloc_Genes_Version_Comparison_by_Cancer_Counts.pdf - Version comparison by cancer (counts)\n")
cat("5. Coloc_Genes_Version_Comparison_by_Cancer_Proportions.pdf - Version comparison by cancer (proportions)\n")

cat("\nGWAS LOCI PLOTS:\n")
cat("6. Coloc_GWAS_Loci_Method_Version_Comparison_Main.pdf - Combined summary plot\n")
cat("7. Coloc_GWAS_Loci_Method_Comparison_by_Cancer_Counts.pdf - Method comparison by cancer (counts)\n")
cat("8. Coloc_GWAS_Loci_Method_Comparison_by_Cancer_Proportions.pdf - Method comparison by cancer (proportions)\n")
cat("9. Coloc_GWAS_Loci_Version_Comparison_by_Cancer_Counts.pdf - Version comparison by cancer (counts)\n")
cat("10. Coloc_GWAS_Loci_Version_Comparison_by_Cancer_Proportions.pdf - Version comparison by cancer (proportions)\n")

cat("\nTABLES:\n")
cat("1. Table_Coloc_Genes_Method_Overlap_by_Cancer.csv - Colocalized genes by method overlap\n")
cat("2. Table_Coloc_Genes_Version_Overlap_by_Cancer.csv - Colocalized genes by version overlap\n")
cat("3. Table_Coloc_GWAS_Loci_Method_Overlap_by_Cancer.csv - GWAS loci by method overlap\n")
cat("4. Table_Coloc_GWAS_Loci_Version_Overlap_by_Cancer.csv - GWAS loci by version overlap\n")

cat("\nAll files saved to:", output_dir, "\n")
