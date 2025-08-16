setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/4_twas/')
require(data.table)
twas_res = fread('all_TWAS_results.txt')
twas_res$TWAS_P = 2*pnorm(-1*abs(twas_res$TWAS_Z))

### Load libraries
library(dplyr)
library(ggplot2)
library(grid)
library(tidyr)
library(paletteer)
library(scales)
library(gridExtra)
library(hexbin)  # For hexagonal binning

twas_res <- twas_res %>%
  mutate(Gene_Tissue = paste(Gene, Tissue, sep = "_"))


# -------------------
# 1. SALMON vs STAR scatter plot
# -------------------
# For method comparison, we need to also group by Version to avoid duplicates
both_methods <- twas_res %>%
  group_by(Gene_Tissue, Version, Method) %>%
  summarise(N_SNPs = mean(N_SNPs, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene_Tissue, Version) %>%
  filter(n_distinct(Method) == 2) %>%
  ungroup() %>%
  pivot_wider(names_from = Method, values_from = N_SNPs, names_prefix = "Method_")

p1 <- ggplot(both_methods, aes(x = Method_SALMON, y = Method_STAR)) +
  geom_hex(bins = 50, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.2) +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "#SNPs in model (Salmon)",
    y = "#SNPs in model (STAR)"
  ) 
# -------------------
# 2-4. Version comparison scatter plots
# -------------------
# For version comparison, we need to also group by Method to avoid duplicates
all_versions <- twas_res %>%
  group_by(Gene_Tissue, Method, Version) %>%
  summarise(N_SNPs = mean(N_SNPs, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene_Tissue, Method) %>%
  filter(n_distinct(Version) == 3) %>%
  ungroup() %>%
  pivot_wider(names_from = Version, values_from = N_SNPs, names_prefix = "v")

# v27 vs v38
p2 <- ggplot(all_versions, aes(x = vv27, y = vv38)) +
  geom_hex(bins = 50, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.2) +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "#SNPs in model (v27)",
    y = "#SNPs in model (v38)"
  ) 

# v27 vs v45
p3 <- ggplot(all_versions, aes(x = vv27, y = vv45)) +
  geom_hex(bins = 50, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.2) +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "#SNPs in model (v27)",
    y = "#SNPs in model (v45)"
  ) 

# v38 vs v45
p4 <- ggplot(all_versions, aes(x = vv38, y = vv45)) +
  geom_hex(bins = 50, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.2) +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "#SNPs in model (v38)",
    y = "#SNPs in model (v45)"
  ) 

# -------------------
# Save plots
# -------------------

require(cowplot)
combined_plot <- plot_grid(p1, p2, p3, p4, ncol = 2,
                           labels = c('A','B','C','D'),
                           label_size = 11)
ggsave("/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/N_SNPs_all_scatter_plots.pdf", 
       combined_plot, width = 8, height = 6)



setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/4_twas/')
require(data.table)
twas_res = fread('all_TWAS_results.txt')
twas_res$TWAS_P = 2*pnorm(-1*abs(twas_res$TWAS_Z))

twas_res = twas_res %>%
  group_by(Gene,Cancer,Version,Method) %>%
  summarize(Gene = first(Gene),
            Cancer = first(Cancer),
            Version = first(Version),
            Method = first(Method),
            Top_GWAS_SNP = first(Top_GWAS_SNP),
            Top_GWAS_P = first(Top_GWAS_P),
            P = ACAT::ACAT(TWAS_P))

twas_res = subset(twas_res,P < 2.5e-6)
twas_res = subset(twas_res,Top_GWAS_P < 5e-8)


# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(paletteer)
library(scales)
library(cowplot)

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

# Prepare TWAS data (assuming twas_res is already processed as shown)
# twas_res should be filtered to P < 2.5e-6 and aggregated by Gene,Cancer,Version,Method
twas_res$Cancer <- cancer_names[twas_res$Cancer]

# Create method overlap analysis
method_overlap_twas <- twas_res %>%
  group_by(Cancer, Version, Gene) %>%
  summarise(
    methods = paste(sort(unique(Method)), collapse = ","),
    n_methods = n_distinct(Method),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_group = case_when(
      methods == "SALMON" ~ "Salmon only",
      methods == "STAR" ~ "STAR only", 
      methods == "SALMON,STAR" ~ "Both",
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

# Create version overlap analysis  
version_overlap_twas <- twas_res %>%
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

# Prepare data for plotting
method_plot_data_twas <- method_overlap_twas %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", gsub("v", "", Version)),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

version_plot_data_twas <- version_overlap_twas %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = toupper(Method),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

# Summary plots - aggregate across cancers
# Method summary data
method_summary_data_twas <- method_plot_data_twas %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Version summary data  
version_summary_data_twas <- version_plot_data_twas %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Plot 1: Method comparison by cancer - stacked counts
p1_cancer <- ggplot(method_plot_data_twas, aes(x = Cancer, y = count, fill = overlap_group)) +
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
    y = "Number of Genes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison by cancer - proportions
p2_cancer <- ggplot(method_plot_data_twas, aes(x = Cancer, y = count, fill = overlap_group)) +
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
    y = "Proportion of Genes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison by cancer - stacked counts
p3_cancer <- ggplot(version_plot_data_twas, aes(x = Cancer, y = count, fill = overlap_group)) +
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
    y = "Number of Genes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison by cancer - proportions
p4_cancer <- ggplot(version_plot_data_twas, aes(x = Cancer, y = count, fill = overlap_group)) +
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
    y = "Proportion of Genes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots (main figures)
# Plot 1 Summary: Method comparison summary - proportions
p1_summary_twas <- ggplot(method_summary_data_twas, aes(x = Version, y = total_count, fill = overlap_group)) +
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
    y = "Proportion of Genes",
    fill = "Method"
  )

# Plot 2 Summary: Version comparison summary - proportions
p2_summary_twas <- ggplot(version_summary_data_twas, aes(x = Method, y = total_count, fill = overlap_group)) +
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
    y = "Proportion of Genes", 
    fill = "GENCODE"
  )

# Create combined main figure
p_twas_main <- plot_grid(p1_summary_twas, p2_summary_twas, ncol = 2)

# Display plots
cat("=== CANCER-BY-CANCER PLOTS (Supplemental) ===\n")
print(p1_cancer)
print(p2_cancer)
print(p3_cancer) 
print(p4_cancer)

cat("\n=== MAIN PLOTS (Summary) ===\n")
print(p1_summary_twas)
print(p2_summary_twas)
print(p_twas_main)

# Save plots
ggsave(plot = p_twas_main,
       filename = 
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/TWAS_Method_Version_Comparison_Main.pdf',
       height = 4,
       width = 7)

ggsave(plot = p1_cancer,
       filename = 
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/TWAS_Method_Comparison_by_Cancer_Counts.pdf',
       height = 8,
       width = 8)

ggsave(plot = p2_cancer,
       filename = 
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/TWAS_Method_Comparison_by_Cancer_Proportions.pdf',
       height = 8,
       width = 8)

ggsave(plot = p3_cancer,
       filename = 
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/TWAS_Version_Comparison_by_Cancer_Counts.pdf',
       height = 8,
       width = 8)

ggsave(plot = p4_cancer,
       filename = 
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/TWAS_Version_Comparison_by_Cancer_Proportions.pdf',
       height = 8,
       width = 8)







# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(paletteer)
library(scales)
library(cowplot)
library(biomaRt)

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

# Get SNP locations using biomaRt
cat("Connecting to biomaRt and retrieving SNP locations...\n")
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Get unique SNPs from twas_res
unique_snps <- unique(twas_res$Top_GWAS_SNP)
cat(paste("Found", length(unique_snps), "unique SNPs to query\n"))

# Query biomaRt for SNP locations (in batches to avoid timeout)
batch_size <- 500
snp_locations <- data.frame()

for(i in seq(1, length(unique_snps), batch_size)) {
  end_idx <- min(i + batch_size - 1, length(unique_snps))
  batch_snps <- unique_snps[i:end_idx]
  
  cat(paste("Processing batch", ceiling(i/batch_size), "of", ceiling(length(unique_snps)/batch_size), "\n"))
  
  batch_result <- getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
    filters = "snp_filter",
    values = batch_snps,
    mart = mart
  )
  
  snp_locations <- rbind(snp_locations, batch_result)
}

# Clean up SNP location data
snp_locations <- snp_locations %>%
  filter(chr_name %in% 1:22) %>%  # Keep only autosomal chromosomes
  mutate(
    chr_name = as.numeric(chr_name),
    position = chrom_start  # Use start position
  ) %>%
  dplyr::select(refsnp_id, chr_name, position) %>%
  distinct()

cat(paste("Retrieved locations for", nrow(snp_locations), "SNPs\n"))

# Merge SNP locations with TWAS results
twas_res_with_loc <- twas_res %>%
  left_join(snp_locations, by = c("Top_GWAS_SNP" = "refsnp_id")) %>%
  filter(!is.na(chr_name) & !is.na(position))  # Remove SNPs without location data

cat(paste("TWAS results with locations:", nrow(twas_res_with_loc), "of", nrow(twas_res), "original rows\n"))

# Define GWAS loci (SNPs within 1MB of each other) grouped by cancer
twas_res_with_loci <- twas_res_with_loc %>%
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

# Create method overlap analysis for GWAS loci
method_overlap_loci <- twas_res_with_loci %>%
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
      methods == "SALMON" ~ "Salmon only",
      methods == "STAR" ~ "STAR only", 
      methods == "SALMON,STAR" ~ "Both",
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
version_overlap_loci <- twas_res_with_loci %>%
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

# Prepare data for plotting
method_plot_data_loci <- method_overlap_loci %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", gsub("v", "", Version)),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

version_plot_data_loci <- version_overlap_loci %>%
  mutate(
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = toupper(Method),
    Cancer_clean = cancer_names[Cancer],
    Cancer_clean = factor(Cancer_clean, levels = cancer_names)
  )

# Summary plots - aggregate across cancers
# Method summary data
method_summary_data_loci <- method_plot_data_loci %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

# Version summary data  
version_summary_data_loci <- version_plot_data_loci %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  )

# Plot 1: Method comparison by cancer - stacked counts
p1_cancer_loci <- ggplot(method_plot_data_loci, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
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
    y = "Number of GWAS Loci\ntagged by TWAS",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison by cancer - proportions
p2_cancer_loci <- ggplot(method_plot_data_loci, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
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
    y = "Proportion of GWAS Loci\ntagged by TWAS",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison by cancer - stacked counts
p3_cancer_loci <- ggplot(version_plot_data_loci, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
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
    y = "Number of GWAS Loci\ntagged by TWAS",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison by cancer - proportions
p4_cancer_loci <- ggplot(version_plot_data_loci, aes(x = Cancer_clean, y = count, fill = overlap_group)) +
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
    y = "Proportion of GWAS Loci\ntagged by TWAS",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots (main figures)
# Plot 1 Summary: Method comparison summary - proportions
p1_summary_loci <- ggplot(method_summary_data_loci, aes(x = Version, y = total_count, fill = overlap_group)) +
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
    y = "Proportion of GWAS Loci\ntagged by TWAS",
    fill = "Method"
  )

# Plot 2 Summary: Version comparison summary - proportions
p2_summary_loci <- ggplot(version_summary_data_loci, aes(x = Method, y = total_count, fill = overlap_group)) +
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
    y = "Proportion of GWAS Loci\ntagged by TWAS", 
    fill = "GENCODE"
  )

# Create combined main figure
p_loci_main <- plot_grid(p1_summary_loci, p2_summary_loci, ncol = 2)

# Print summary statistics
cat("\n=== GWAS LOCI SUMMARY ===\n")
total_loci_by_cancer <- twas_res_with_loci %>%
  group_by(Cancer) %>%
  summarise(unique_loci = n_distinct(gwas_locus), .groups = "drop") %>%
  mutate(Cancer_clean = cancer_names[Cancer])

print(total_loci_by_cancer)

cat(paste("\nTotal unique GWAS loci across all cancers:", 
          n_distinct(twas_res_with_loci$gwas_locus), "\n"))

# Save plots to specified directory
output_dir <- "/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/"

ggsave(plot = p_loci_main,
       filename = paste0(output_dir, "TWAS_GWAS_Loci_Method_Version_Comparison_Main.pdf"),
       height = 4,
       width = 7)

ggsave(plot = p1_cancer_loci,
       filename = paste0(output_dir, "TWAS_GWAS_Loci_Method_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p2_cancer_loci,
       filename = paste0(output_dir, "TWAS_GWAS_Loci_Method_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p3_cancer_loci,
       filename = paste0(output_dir, "TWAS_GWAS_Loci_Version_Comparison_by_Cancer_Counts.pdf"),
       height = 8,
       width = 8)

ggsave(plot = p4_cancer_loci,
       filename = paste0(output_dir, "TWAS_GWAS_Loci_Version_Comparison_by_Cancer_Proportions.pdf"),
       height = 8,
       width = 8)

cat("\nPlots saved to:", output_dir, "\n")



# Export tables to CSV files using the correct gene analysis data
output_dir <- "/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/"

# Table 1: Number of TWAS genes per cancer - Method overlap (uses method_plot_data_twas)
twas_genes_method_table <- method_plot_data_twas %>%
  dplyr::select(Cancer_clean, Version, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean,
    `GENCODE Version` = Version
  ) %>%
  arrange(`Cancer Type`, `GENCODE Version`)

# Table 2: Number of TWAS genes per cancer - Version overlap (uses version_plot_data_twas)  
twas_genes_version_table <- version_plot_data_twas %>%
  dplyr::select(Cancer_clean, Method, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean
  ) %>%
  arrange(`Cancer Type`, Method)

# Table 3: Number of GWAS loci per cancer - Method overlap (uses method_plot_data_loci)
gwas_loci_method_table <- method_plot_data_loci %>%
  dplyr::select(Cancer_clean, Version, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean,
    `GENCODE Version` = Version
  ) %>%
  arrange(`Cancer Type`, `GENCODE Version`)

# Table 4: Number of GWAS loci per cancer - Version overlap (uses version_plot_data_loci)
gwas_loci_version_table <- version_plot_data_loci %>%
  dplyr::select(Cancer_clean, Method, overlap_group, count) %>%
  pivot_wider(names_from = overlap_group, values_from = count, values_fill = 0) %>%
  rename(
    `Cancer Type` = Cancer_clean
  ) %>%
  arrange(`Cancer Type`, Method)

# Write tables to CSV files
write.csv(twas_genes_method_table, 
          paste0(output_dir, "Table_TWAS_Genes_Method_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(twas_genes_version_table,
          paste0(output_dir, "Table_TWAS_Genes_Version_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(gwas_loci_method_table,
          paste0(output_dir, "Table_GWAS_Loci_Method_Overlap_by_Cancer.csv"),
          row.names = FALSE)

write.csv(gwas_loci_version_table,
          paste0(output_dir, "Table_GWAS_Loci_Version_Overlap_by_Cancer.csv"),
          row.names = FALSE)





# Calculate genes per GWAS locus that are NOT picked up by both methods
genes_per_locus_method_discordant <- twas_res_with_loci %>%
  group_by(Cancer, gwas_locus, Gene) %>%
  summarise(
    methods = paste(sort(unique(Method)), collapse = ","),
    n_methods = n_distinct(Method),
    .groups = "drop"
  ) %>%
  # Keep only genes that are method-specific (not in both)
  filter(n_methods == 1) %>%
  group_by(Cancer, gwas_locus) %>%
  summarise(
    discordant_genes_per_locus = n(),
    .groups = "drop"
  ) %>%
  mutate(Cancer_clean = cancer_names[Cancer])

# Calculate genes per GWAS locus that are NOT in all three versions
genes_per_locus_version_discordant <- twas_res_with_loci %>%
  group_by(Cancer, gwas_locus, Gene) %>%
  summarise(
    versions = paste(sort(unique(Version)), collapse = ","),
    n_versions = n_distinct(Version),
    .groups = "drop"
  ) %>%
  # Keep only genes that are not in all three versions
  filter(n_versions < 3) %>%
  group_by(Cancer, gwas_locus) %>%
  summarise(
    discordant_genes_per_locus = n(),
    .groups = "drop"
  ) %>%
  mutate(Cancer_clean = cancer_names[Cancer])

# Box plot 1: Method discordant genes per locus (one box per cancer)
p_method_boxplot <- ggplot(genes_per_locus_method_discordant, aes(x = Cancer_clean, y = discordant_genes_per_locus)) +
  stat_summary(
    fun.data = function(y) {
      qs <- quantile(y, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      data.frame(
        middle = qs[3],                         # median (50%)
        ymin = qs[1], ymax = qs[5],            # 10–90%
        lower = qs[2], upper = qs[4]           # 25–75%
      )
    },
    geom = "boxplot", alpha = 0.7, fill = "#E31A1C"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Method-discordant genes per GWAS locus"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# Box plot 2: Version discordant genes per locus (one box per cancer)
p_version_boxplot <- ggplot(genes_per_locus_version_discordant, aes(x = Cancer_clean, y = discordant_genes_per_locus)) +
  stat_summary(
    fun.data = function(y) {
      qs <- quantile(y, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      data.frame(
        middle = qs[3],                         # median (50%)
        ymin = qs[1], ymax = qs[5],            # 10–90%
        lower = qs[2], upper = qs[4]           # 25–75%
      )
    },
    geom = "boxplot", alpha = 0.7, fill = "#1B9E77"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  ) +
  labs(
    x = "Cancer Type",
    y = "Version-discordant genes per GWAS locus"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# Print summary statistics
cat("\n=== DISCORDANT GENES PER GWAS LOCUS SUMMARY ===\n")
cat("Method-discordant genes summary:\n")
method_summary <- genes_per_locus_method_discordant %>%
  group_by(Cancer_clean) %>%
  summarise(
    median_genes = median(discordant_genes_per_locus),
    mean_genes = round(mean(discordant_genes_per_locus), 2),
    max_genes = max(discordant_genes_per_locus),
    n_loci = n(),
    .groups = "drop"
  )
print(method_summary)

cat("\nVersion-discordant genes summary:\n")
version_summary <- genes_per_locus_version_discordant %>%
  group_by(Cancer_clean) %>%
  summarise(
    median_genes = median(discordant_genes_per_locus),
    mean_genes = round(mean(discordant_genes_per_locus), 2),
    max_genes = max(discordant_genes_per_locus),
    n_loci = n(),
    .groups = "drop"
  )
print(version_summary)

# Combine both discordance measures into one dataset
combined_discordance <- bind_rows(
  genes_per_locus_method_discordant %>%
    dplyr::select(Cancer_clean, discordant_genes_per_locus) %>%
    mutate(Discordance_Type = "Method"),
  genes_per_locus_version_discordant %>%
    dplyr::select(Cancer_clean, discordant_genes_per_locus) %>%
    mutate(Discordance_Type = "Annotation")
) %>%
  mutate(
    Discordance_Type = factor(Discordance_Type, levels = c("Method", "Annotation"))
  )

# Combined box plot: Both discordance types (two boxes per cancer)
p_combined_boxplot <- ggplot(combined_discordance, aes(x = Cancer_clean, y = discordant_genes_per_locus, fill = Discordance_Type)) +
  stat_summary(
    fun.data = function(y) {
      qs <- quantile(y, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      data.frame(
        middle = qs[3],                         # median (50%)
        ymin = qs[1], ymax = qs[5],            # 10–90%
        lower = qs[2], upper = qs[4]           # 25–75%
      )
    },
    geom = "boxplot", alpha = 0.7, position = position_dodge(width = 0.6),width = .5
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(2, "mm"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "Cancer Type",
    y = "Discordant genes per GWAS locus",
    fill = "Discordance Type"
  ) +
  scale_fill_manual(values = c("Method" = "#E31A1C", "Annotation" = "#1B9E77")) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# Print summary statistics
cat("\n=== COMBINED DISCORDANT GENES PER GWAS LOCUS SUMMARY ===\n")
combined_summary <- combined_discordance %>%
  group_by(Cancer_clean, Discordance_Type) %>%
  summarise(
    median_genes = median(discordant_genes_per_locus),
    mean_genes = round(mean(discordant_genes_per_locus), 2),
    max_genes = max(discordant_genes_per_locus),
    n_loci = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Discordance_Type, values_from = c(median_genes, mean_genes, max_genes, n_loci),
              names_sep = "_")
print(combined_summary)

# Display plot
print(p_combined_boxplot)

# Save combined plot
output_dir <- "/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/"

ggsave(plot = p_combined_boxplot,
       filename = paste0(output_dir, "Combined_Discordant_Genes_per_GWAS_Locus.pdf"),
       height = 6,
       width = 6)

# Also save individual plots for reference
ggsave(plot = p_method_boxplot,
       filename = paste0(output_dir, "Method_Discordant_Genes_per_GWAS_Locus.pdf"),
       height = 6,
       width = 6)

ggsave(plot = p_version_boxplot,
       filename = paste0(output_dir, "Version_Discordant_Genes_per_GWAS_Locus.pdf"),
       height = 6,
       width = 6)

cat("\nBox plots saved to:", output_dir, "\n")