library(jsonlite)
library(dplyr)
library(tidyr)

library(jsonlite)

# --- Read sequencing.json ---
txt <- readLines("/rsrch5/home/epi/stbresnahan/scratch/GTEx/sequencing.json",
                 warn = FALSE, encoding = "UTF-8")

# Clean invalid UTF8
txt_clean <- iconv(txt, from = "", to = "UTF-8", sub = "")

seq_json <- fromJSON(txt_clean, flatten = TRUE)
df <- as.data.frame(seq_json[["data"]])

# Filter for STAR-aligned BAMs (your sequencing criteria)
bam_files <- df$file_name[
  df$data_type == "Aligned Reads" &
    df$file_format == "bam" &
    df$data_category == "Sequencing Reads"
]

df$file_name[
  df$data_type == "Unaligned Reads" &
    df$file_format == "bam" &
    df$data_category == "Sequencing Reads"
]

bam_df <- data.frame(
  file_name = bam_files,
  sample_id = sub("\\..*$", "", bam_files),  # everything before the first period
  stringsAsFactors = FALSE
)

# --- Read sample attributes ---
sample_attributes <- read.delim(
  "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx_v8_sample_attributes.txt",
  stringsAsFactors = FALSE
)

# Only keep SAMPID and tissue columns
sample_tissue <- sample_attributes %>%
  select(SAMPID, tissue)

bam_df <- bam_df[bam_df$sample_id%in%sample_tissue$SAMPID,]

bam_df <- bam_df %>%
  left_join(sample_tissue, by = c("sample_id" = "SAMPID"))


# --- Read TSV list created by Bash script ---
bam_list <- read.delim("/rsrch5/home/epi/stbresnahan/scratch/GTEx/GTE_bam_list.tsv",
                       stringsAsFactors = FALSE)

# --- Check if each bam_file from the TSV exists in sequencing.json ---
bam_list$match_in_seq_json <- bam_list$bam_file %in% bam_df$file_name

# --- Result ---
head(bam_list)

# --- Merge tissue into bam_list ---
bam_list <- bam_list %>%
  left_join(sample_tissue, by = c("sample_id" = "SAMPID"))

# Reorder columns if you like
bam_list <- bam_list %>%
  select(tissue, everything())

summary_table <- bam_list %>%
  group_by(tissue) %>%
  summarise(
    total_samples = n(),
    true_matches = sum(match_in_seq_json, na.rm = TRUE),
    false_matches = sum(!match_in_seq_json, na.rm = TRUE)
  ) %>%
  arrange(tissue)

summary_table
