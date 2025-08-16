setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/4_twas/')

res = readRDS('twas_counts.RDS')
res$Gene = sapply(strsplit(res$File, '/'), function(x) x[10])
res$Gene = sapply(strsplit(res$Gene, '_TWAS'), function(x) x[1])

library(ggplot2)
library(dplyr)
library(tidyr)
library(paletteer)
library(scales)
library(purrr)

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

# Read sample sizes and create tissue order
n_gtex = data.table::fread('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/3_results/N_GTEx.tsv')
n_gtex = n_gtex[order(n_gtex$N, decreasing = T),]
tissue_order <- tissue_abbrev[n_gtex$Tissue]

# Method comparison: Create overlap groups (STAR only, Salmon only, Both)
method_overlap <- res %>%
  group_by(Tissue, Version, Gene) %>%
  summarise(
    methods = list(Method),
    .groups = "drop"
  ) %>%
  mutate(
    has_star = map_lgl(methods, ~ "STAR" %in% .x),
    has_salmon = map_lgl(methods, ~ "Salmon" %in% .x),
    overlap_group = case_when(
      has_star & has_salmon ~ "Both",
      has_star & !has_salmon ~ "STAR only",
      !has_star & has_salmon ~ "Salmon only"
    )
  ) %>%
  group_by(Tissue, Version, overlap_group) %>%
  summarise(count = n(), .groups = "drop")

# Keep version information for method comparison plots
method_plot_data <- method_overlap %>%
  mutate(
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("STAR only", "Salmon only", "Both")),
    Version = paste0(Version)
  )

# Version comparison: Create overlap groups 
version_overlap <- res %>%
  group_by(Tissue, Method, Gene) %>%
  summarise(
    versions = list(as.character(Version)),
    .groups = "drop"
  ) %>%
  mutate(
    has_v27 = map_lgl(versions, ~ "v27" %in% .x),
    has_v38 = map_lgl(versions, ~ "v38" %in% .x),
    has_v45 = map_lgl(versions, ~ "v45" %in% .x),
    overlap_group = case_when(
      has_v27 & has_v38 & has_v45 ~ "All three",
      has_v27 & has_v38 & !has_v45 ~ "v27 & v38",
      has_v27 & !has_v38 & has_v45 ~ "v27 & v45", 
      !has_v27 & has_v38 & has_v45 ~ "v38 & v45",
      has_v27 & !has_v38 & !has_v45 ~ "v27 only",
      !has_v27 & has_v38 & !has_v45 ~ "v38 only",
      !has_v27 & !has_v38 & has_v45 ~ "v45 only"
    )
  ) %>%
  group_by(Tissue, Method, overlap_group) %>%
  summarise(count = n(), .groups = "drop")

# Keep method information for version comparison plots
version_plot_data <- version_overlap %>%
  mutate(
    tissue_abbrev = tissue_abbrev[Tissue],
    tissue_abbrev = factor(tissue_abbrev, levels = tissue_order),
    overlap_group = factor(overlap_group, levels = c("v27 only", "v38 only", "v45 only", 
                                                     "v27 & v38", "v27 & v45", "v38 & v45", "All three")),
    Method = toupper(Method)
  )

# Plot 1: Method comparison - Total TWAS genes by overlap group, faceted by version
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
    y = expression("Number of genes with"~R^2 >= 0.01),
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C"))

# Plot 2: Method comparison - Relative proportion within each version
p2 <- ggplot(method_plot_data, aes(x = tissue_abbrev, y = count, fill = overlap_group)) +
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
    y = expression("Proportion of genes with"~R^2 >= 0.01),
    fill = "Method"
  ) +
  scale_fill_manual(values = c("STAR only" = "#E31A1C", "Salmon only" = "#1F78B4", "Both" = "#33A02C")) +
  scale_y_continuous(labels = scales::percent_format())

# Plot 3: Version comparison - Total TWAS genes by overlap group, faceted by method
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
    y = expression("Number of genes with"~R^2 >= 0.01),
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7))

# Plot 4: Version comparison - Relative proportion within each method
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
    y = expression("Proportion of genes with"~R^2 >= 0.01),
    fill = "GENCODE"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_jama", n = 7)) +
  scale_y_continuous(labels = scales::percent_format())

# Summary plots - aggregate across tissues
# Method summary data
method_summary_data <- method_plot_data %>%
  group_by(Version, overlap_group) %>%
  summarise(
    total_count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Version summary data  
version_summary_data <- version_plot_data %>%
  group_by(Method, overlap_group) %>%
  summarise(
    total_count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(total_count_thousands = total_count / 1000)

# Plot 1 Summary: Method comparison summary - Total TWAS genes
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
    y = expression("Number of genes with"~R^2 >= 0.01~"("~x10^3~")"),
    fill = "Method"
  )

# Plot 2 Summary: Method comparison summary - Relative proportion
p2_summary <- ggplot(method_summary_data, aes(x = Version, y = total_count, fill = overlap_group)) +
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
    y = expression("Proportion of genes with"~R^2 >= 0.01),
    fill = "Method"
  )

# Plot 3 Summary: Version comparison summary - Total TWAS genes
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
    y = expression("Number of genes with"~R^2 >= 0.01~"("~x10^3~")"),
    fill = "GENCODE"
  )

# Plot 4 Summary: Version comparison summary - Relative proportion
p4_summary <- ggplot(version_summary_data, aes(x = Method, y = total_count, fill = overlap_group)) +
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
    y = expression("Proportion of genes with"~R^2 >= 0.01), 
    fill = "GENCODE"
  )

# Method consistency for TWAS genes - proportion of "Both" out of total TWAS genes per tissue
method_consistency_twas <- method_plot_data %>%
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

# Version consistency for TWAS genes - proportion of "All three" out of total TWAS genes per tissue  
version_consistency_twas <- version_plot_data %>%
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
    Method = factor(Method, levels = c("STAR", "SALMON"))
  )

# TWAS Plot 1: Method consistency (Both methods) vs Sample Size, colored by Version
p1_twas <- ggplot(method_consistency_twas, aes(x = N, y = 1-consistency_proportion, color = Version)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid")  +
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
    legend.position = c(0.02, 0.3),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "Gene prediction discordance %\nacross method",
    color = "GENCODE"
  ) +
  scale_y_continuous(labels = scales::percent_format())

# TWAS Plot 2: Version consistency (All three versions) vs Sample Size, colored by Method
p2_twas <- ggplot(version_consistency_twas, aes(x = N, y = 1-consistency_proportion, color = Method)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = F, linetype = "solid") +
  scale_color_manual(values = c("STAR" = "#E31A1C", "SALMON" = "#1F78B4")) +
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
    legend.position = c(0.02, 0.3),
    legend.justification = c(0, 1)
  ) +
  labs(
    x = "N",
    y = "Gene prediction discordance %\nacross annotation",
    color = "Method"
  ) +
  scale_y_continuous(labels = scales::percent_format())

# Display tissue-by-tissue plots
cat("=== TISSUE-BY-TISSUE PLOTS ===\n")
print(p1)
print(p2)
print(p3) 
print(p4)

# Display summary plots
cat("\n=== SUMMARY PLOTS ===\n")
print(p1_summary)
print(p2_summary)
print(p3_summary)
print(p4_summary)

# Display TWAS plots
print(p1_twas)
print(p2_twas)

require(cowplot)
p_tot = plot_grid(p2_summary,
                  p4_summary,
                  ncol = 2)
ggsave(plot = p_tot,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Bottom.pdf',
       height = 4,
       width = 7)

ggsave(plot = p1_twas,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Scatter_TWASPred_Method.pdf',
       height = 4,
       width = 2.5)
ggsave(plot = p2_twas,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/Fig1Scatter_TWASPred_Version.pdf',
       height = 4,
       width = 2.5)

ggsave(plot = p1,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommonTWASPred_Method_Tissue_N.pdf',
       height = 8,
       width = 8)
ggsave(plot = p2,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommonTWASPred_Method_Tissue_Proportion.pdf',
       height = 8,
       width = 8)
ggsave(plot = p3,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommonTWASPred_Version_Tissue_N.pdf',
       height = 8,
       width = 8)
ggsave(plot = p4,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/Figures/SuppFig1_CommonTWASPred_Version_Tissue_Proportion.pdf',
       height = 8,
       width = 8)

data.table::fwrite(merge(method_overlap,n_gtex,by='Tissue'),
                   '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/SupplementalTable5.csv',
                   col.names=T,
                   row.names=F,
                   sep=',',
                   quote=F)
data.table::fwrite(merge(version_overlap,n_gtex,by='Tissue'),
                   '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Summer 2025/GTEx GENCODE Comp/Manuscript/SupplementalTable6.csv',
                   col.names=T,
                   row.names=F,
                   sep=',',
                   quote=F)
