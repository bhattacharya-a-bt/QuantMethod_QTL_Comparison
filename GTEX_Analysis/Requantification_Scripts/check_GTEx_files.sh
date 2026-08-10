#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before running.
#   sample_IDs  GTEx v8 sample attributes file (tab-delimited; column 1 is the
#               sample ID)
#   bamdir      directory holding the source GTEx BAM files
#   out_file    output sample_id -> bam_file lookup table
# ============================================================================
sample_IDs="${sample_IDs:-/path/to/GTEx_v8/GTEx_v8_sample_attributes.txt}"
bamdir="${bamdir:-/path/to/GTEx/SourceFiles/Bam}"
out_file="${out_file:-/path/to/scratch/GTEx/GTE_bam_list.tsv}"

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
