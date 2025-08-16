setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/3_results/')
res = readRDS('eGeneCounts.RDS')

n_gtex = data.table::fread('N_GTEx.tsv')
n_gtex = n_gtex[order(n_gtex$N,decreasing = T),]

library(ggplot2)
library(dplyr)
library(tidyr)
library(paletteer)
library(scales)

# Create tissue abbreviation mapping
tissue_abbrev <- c(
  "Adipose_Subcutaneous" = "AdipSub",
  "Adipose_Visceral_Omentum" = "AdipVisc", 
  "Adrenal_Gland" = "Adrenal",
  "Artery_Aorta" = "ArtAorta",
  "Artery_Coronary" = "ArtCor",
  "Artery_Tibial" = "ArtTib",
  "Brain_Amygdala" = "BrAmyg",
  "Brain_Anterior_cingulate_cortex_BA24" = "BrACC",
  "Brain_Caudate_basal_ganglia" = "BrCaud",
  "Brain_Cerebellar_Hemisphere" = "BrCereb",
  "Brain_Cerebellum" = "BrCerebel",
  "Brain_Cortex" = "BrCortex",
  "Brain_Frontal_Cortex_BA9" = "BrFront",
  "Brain_Hippocampus" = "BrHippo",
  "Brain_Hypothalamus" = "BrHypo",
  "Brain_Nucleus_accumbens_basal_ganglia" = "BrNucAcc",
  "Brain_Putamen_basal_ganglia" = "BrPut",
  "Brain_Spinal_cord_cervical_c-1" = "BrSpinal",
  "Brain_Substantia_nigra" = "BrSubNig",
  "Breast_Mammary_Tissue" = "Breast",
  "Cells_Cultured_fibroblasts" = "Fibroblast",
  "Cells_EBV-transformed_lymphocytes" = "Lymph",
  "Colon_Sigmoid" = "ColSig",
  "Colon_Transverse" = "ColTrans",
  "Esophagus_Gastroesophageal_Junction" = "EsoGEJ",
  "Esophagus_Mucosa" = "EsoMuc",
  "Esophagus_Muscularis" = "EsoMus",
  "Heart_Atrial_Appendage" = "HeartAtr",
  "Heart_Left_Ventricle" = "HeartLV",
  "Liver" = "Liver",
  "Lung" = "Lung",
  "Minor_Salivary_Gland" = "Salivary",
  "Muscle_Skeletal" = "Muscle",
  "Nerve_Tibial" = "Nerve",
  "Ovary" = "Ovary",
  "Pancreas" = "Pancreas",
  "Pituitary" = "Pituitary",
  "Prostate" = "Prostate",
  "Skin_Not_Sun_Exposed_Suprapubic" = "SkinNoSun",
  "Skin_Sun_Exposed_Lower_leg" = "SkinSun",
  "Small_Intestine_Terminal_Ileum" = "SmIntestine",
  "Spleen" = "Spleen",
  "Stomach" = "Stomach",
  "Testis" = "Testis",
  "Thyroid" = "Thyroid",
  "Uterus" = "Uterus",
  "Vagina" = "Vagina",
  "Whole_Blood" = "Blood"
)

# Create tissue order based on n_gtex (largest to smallest sample size)
tissue_order <- tissue_abbrev[n_gtex$Tissue]

# Method comparison: Create overlap groups (STAR only, Salmon only, Both)
method_overlap <- res$method_overlap

# Keep version information for method comparison plots
method_plot_data <- method_overlap %>%
  mutate(
    proportion = count / total_genes,
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", Version)
  )

# Version comparison: Create overlap groups 
version_overlap <- res$version_overlap

# Keep method information for version comparison plots
version_plot_data <- version_overlap %>%
  mutate(
    proportion = count / total_genes,
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = tools::toTitleCase(as.character(Method))
  )

# Plot 1: Method comparison - Total eGenes by overlap group, faceted by version
p1 <- ggplot(method_plot_data, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Number of eGenes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison - Proportion of total genes, faceted by version
p2 <- ggplot(method_plot_data, aes(x = tissue_abbrev, y = proportion, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Proportion of Total Genes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

version_plot_data$Method = ifelse(version_plot_data$Method == 'Star',
                                  'STAR',
                                  'Salmon')
# Plot 3: Version comparison - Total eGenes by overlap group, faceted by method
p3 <- ggplot(version_plot_data, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Number of eGenes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison - Proportion of total genes, faceted by method
p4 <- ggplot(version_plot_data, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Proportion of Total Genes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Summary plots - aggregate across tissues
# Method summary data
method_summary_data <- method_plot_data %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Version summary data  
version_summary_data <- version_plot_data %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Plot 1 Summary: Method comparison summary - Total eGenes
p1_summary <- ggplot(method_summary_data, aes(x = Version, y = total_count_thousands, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
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
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(breaks = pretty_breaks(), labels = comma) +
  labs(
    x = "GENCODE",
    y = "Total Number of eGenes (in thousands)",
    fill = "Method"
  )

# Plot 2 Summary: Method comparison summary - Average proportion
p2_summary <- ggplot(method_summary_data, aes(x = Version, y = avg_proportion, fill = overlap_group)) +
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
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  labs(
    x = "GENCODE", 
    y = "Average Proportion of Total Genes",
    fill = "Method"
  )

# Plot 3 Summary: Version comparison summary - Total eGenes
p3_summary <- ggplot(version_summary_data, aes(x = Method, y = total_count_thousands, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
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
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(breaks = pretty_breaks(), labels = comma) +
  labs(
    x = "Method",
    y = "Total Number of eGenes (in thousands)",
    fill = "GENCODE"
  )

# Plot 4 Summary: Version comparison summary - Average proportion  
p4_summary <- ggplot(version_summary_data, aes(x = Method, y = avg_proportion, fill = overlap_group)) +
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
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  labs(
    x = "Method",
    y = "Average Proportion of Total Genes", 
    fill = "GENCODE"
  )

# COMMON DATA PROCESSING AND PLOTS
# Method comparison: Create overlap groups (STAR only, Salmon only, Both) for COMMON
method_overlap_common <- res$method_overlap_common

# Keep version information for method comparison plots - COMMON
method_plot_data_common <- method_overlap_common %>%
  mutate(
    proportion = count / total_genes,
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0("v", Version)
  )

# Version comparison: Create overlap groups for COMMON
version_overlap_common <- res$version_overlap_common

# Keep method information for version comparison plots - COMMON
version_plot_data_common <- version_overlap_common %>%
  mutate(
    proportion = count / total_genes,
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = tools::toTitleCase(as.character(Method))
  )

# COMMON PLOTS - tissue by tissue
# Plot 1: Method comparison - Total eGenes by overlap group, faceted by version - COMMON
p1_common <- ggplot(method_plot_data_common, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Number of eGenes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison - Relative proportion within each version - COMMON
p2_common <- ggplot(method_plot_data_common, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Proportion of eGenes",
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison - Total eGenes by overlap group, faceted by method - COMMON
version_plot_data_common$Method = ifelse(version_plot_data_common$Method == 'Star',
                                         'STAR','Salmon')
p3_common <- ggplot(version_plot_data_common, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Number of eGenes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison - Relative proportion within each method - COMMON
p4_common <- ggplot(version_plot_data_common, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    x = "Tissue",
    y = "Proportion of eGenes",
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots - aggregate across tissues - COMMON
# Method summary data - COMMON
method_summary_data_common <- method_plot_data_common %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Version summary data - COMMON
version_summary_data_common <- version_plot_data_common %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    avg_proportion = mean(proportion),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Plot 1 Summary: Method comparison summary - Total eGenes - COMMON
p1_summary_common <- ggplot(method_summary_data_common, aes(x = Version, y = total_count_thousands, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
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
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(breaks = pretty_breaks(), labels = comma) +
  labs(
    x = "GENCODE",
    y = "Total Number of eGenes (in thousands)",
    fill = "Method"
  )

# Plot 2 Summary: Method comparison summary - Relative proportion - COMMON
p2_summary_common <- ggplot(method_summary_data_common, aes(x = Version, y = total_count, fill = overlap_group)) +
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
    legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0)
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "GENCODE", 
    y = "Proportion of eGenes",
    fill = "Method"
  )

# Plot 3 Summary: Version comparison summary - Total eGenes - COMMON
p3_summary_common <- ggplot(version_summary_data_common, aes(x = Method, y = total_count_thousands, fill = overlap_group)) +
  geom_bar(stat = "identity", position = "stack") +
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
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(breaks = pretty_breaks(), labels = comma) +
  labs(
    x = "Method",
    y = "Total Number of eGenes (in thousands)",
    fill = "GENCODE"
  )

# Plot 4 Summary: Version comparison summary - Relative proportion - COMMON
p4_summary_common <- ggplot(version_summary_data_common, aes(x = Method, y = total_count, fill = overlap_group)) +
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
    legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
    legend.margin=margin(0,0,0,0)
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Method",
    y = "Proportion of eGenes", 
    fill = "GENCODE"
  )

# eGenes scatterplots (from your existing common data)

# Method consistency for eGenes - proportion of "Both" out of total eGenes per tissue
method_consistency_egenes <- method_plot_data_common %>%
  group_by(Tissue, Version) %>%
  summarise(
    both_count = sum(count[overlap_group == "Both"], na.rm = TRUE),
    total_count = sum(count),
    consistency_proportion = both_count / total_count,
    .groups = "drop"
  ) %>%
  left_join(n_gtex, by = "Tissue") %>%
  filter(!is.na(N)) %>%
  mutate(
    tissue_abbrev = tissue_abbrev[Tissue],
    Version = factor(Version, levels = c("v27", "v38", "v45"))
  )

# Version consistency for eGenes - proportion of "All three" out of total eGenes per tissue  
version_consistency_egenes <- version_plot_data_common %>%
  group_by(Tissue, Method) %>%
  summarise(
    all_three_count = sum(count[overlap_group == "All three"], na.rm = TRUE),
    total_count = sum(count),
    consistency_proportion = all_three_count / total_count,
    .groups = "drop"
  ) %>%
  left_join(n_gtex, by = "Tissue") %>%
  filter(!is.na(N)) %>%
  mutate(
    tissue_abbrev = tissue_abbrev[Tissue],
    Method = factor(Method, levels = c("STAR", "Salmon"))
  )

# eGenes Plot 1: Method consistency (Both methods) vs Sample Size, colored by Version
p1_egenes <- ggplot(method_consistency_egenes, aes(x = N, y = 1-consistency_proportion, color = Version)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid") +
  scale_color_manual(values = paletteer_d("ggsci::default_jama", n = 3)) +
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
    legend.position = c(0.5, 0.75),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "eGene discordance %\nacross method",
    color = "GENCODE"
  ) +
  scale_y_continuous(labels = scales::percent_format())

# eGenes Plot 2: Version consistency (All three versions) vs Sample Size, colored by Method
p2_egenes <- ggplot(version_consistency_egenes, aes(x = N, y = 1-consistency_proportion, color = Method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid") +
  scale_color_manual(values = c("STAR" = "#E31A1C", "Salmon" = "#1F78B4")) +
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
    legend.position = c(0.6, 0.8),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "eGene discordance %\nacross annotation",
    color = "Method"
  ) +
  scale_y_continuous(labels = scales::percent_format())

# Display main plots (tissue-by-tissue, supplemental)
cat("=== SUPPLEMENTAL PLOTS (by Tissue) - TOTAL ===\n")
print(p1)
print(p2)
print(p3) 
print(p4)

# Display summary plots (main figures)
cat("\n=== MAIN PLOTS (Summary) - TOTAL ===\n")
print(p1_summary)
print(p2_summary)
print(p3_summary)
print(p4_summary)

# Display COMMON plots
cat("\n=== SUPPLEMENTAL PLOTS (by Tissue) - COMMON ===\n")
print(p1_common)
print(p2_common)
print(p3_common) 
print(p4_common)

cat("\n=== MAIN PLOTS (Summary) - COMMON ===\n")
print(p1_summary_common)
print(p2_summary_common)
print(p3_summary_common)
print(p4_summary_common)

# Display scatter plots
print(p1_egenes)
print(p2_egenes)

require(cowplot)
p_tot = plot_grid(p2_summary_common,
                  p4_summary_common,
                  ncol = 2)
ggsave(plot = p_tot,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Top.pdf',
       height = 4,
       width = 7)

ggsave(plot = p1_egenes,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Scatter_eGene_Method.pdf',
       height = 4,
       width = 2.5)
ggsave(plot = p2_egenes,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Scatter_eGene_Version.pdf',
       height = 4,
       width = 2.5)

ggsave(plot = p1_common,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommoneGene_Method_Tissue_N.pdf',
       height = 8,
       width = 8)
ggsave(plot = p2_common,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommoneGene_Method_Tissue_Proportion.pdf',
       height = 8,
       width = 8)
ggsave(plot = p3_common,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommoneGene_Version_Tissue_N.pdf',
       height = 8,
       width = 8)
ggsave(plot = p4_common,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommoneGene_Version_Tissue_Proportion.pdf',
       height = 8,
       width = 8)

data.table::fwrite(res$method_overlap_common,
                   '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/SupplementalTable3.csv',
                   col.names=T,
                   row.names=F,
                   sep=',',
                   quote=F)
data.table::fwrite(res$version_overlap_common,
                   '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/SupplementalTable4.csv',
                   col.names=T,
                   row.names=F,
                   sep=',',
                   quote=F)
