#!/bin/bash

# --- Paths ---
sample_IDs="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx_v8_sample_attributes.txt"
bamdir="/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam"
out_file="/rsrch5/home/epi/stbresnahan/scratch/GTEx/GTE_bam_list.tsv"

# --- Prepare output file ---
echo -e "sample_id\tbam_file" > "$out_file"

# --- Loop over all sample IDs (ignore tissue) ---
awk 'NR>1 {print $1}' "$sample_IDs" | while IFS= read -r FILE; do
    MATCH=$(ls "${bamdir}/${FILE}"*.bam 2>/dev/null | head -n 1)
    
    if [[ -n "$MATCH" ]]; then
        echo -e "${FILE}\t$(basename "$MATCH")" >> "$out_file"
    fi
done

echo "GTE_bam_list.tsv created at $out_file"
