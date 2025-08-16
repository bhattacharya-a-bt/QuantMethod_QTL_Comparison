#!/bin/bash
#BSUB -J "concat_twas[1-48]"
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/concat_twas_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/concat_twas_%I.err
#BSUB -n 1
#BSUB -M 40
#BSUB -R "rusage[mem=40G]"
#BSUB -W 24:00
#BSUB -q medium

# Output directory and log setup
mkdir -p /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs

# List the 48 subfolders
SUBFOLDERS=(/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare/*)

# Pick the subfolder for this job
SUBFOLDER=${SUBFOLDERS[$((LSB_JOBINDEX-1))]}

echo "Job ${LSB_JOBINDEX}: Processing folder $SUBFOLDER"

# Output for this job
OUT=/rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/tmp_concat_${LSB_JOBINDEX}.txt
: > "$OUT"

# Concatenate all matching files in this subfolder
find "$SUBFOLDER" -type f -name "*.txt.gz" \
  | grep "/twas/" \
  | while read -r f; do
      zcat "$f" | tail -n +2 >> "$OUT"
    done