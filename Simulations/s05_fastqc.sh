#!/bin/bash

####################################################################################
# Usage:
#   bash s05_fastqc.sh <PASS> <PARAM_ROW_READS> <GENO_PASS>
####################################################################################

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before submitting the job.
#   BASE_DIR   base project directory; the pass number is appended internally
#   CONDA_BIN  path to the conda binary
# ============================================================================
BASE_DIR="${BASE_DIR:-/path/to/scratch/GTEx_gencode_comp}"
CONDA_BIN="${CONDA_BIN:-/path/to/miniforge3/bin/conda}"

module load fastqc

eval "$(${CONDA_BIN} shell.bash hook)"
conda activate fastqc-0.11.9

PASS=$1
PARAM_ROW_READS=$2
GENO_PASS=$3
i=${LSB_JOBINDEX}
sample_file="${BASE_DIR}/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids"
sample_name=$(awk -v row=$i 'NR == row {print $1}' "$sample_file")

# path to your paramspace for reads file
psr_file="${BASE_DIR}/pass${PASS}/files_for_analysis/parameter_space_reads.txt"

# extract the paired_end status from the specified row and column
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

reads_dir="${BASE_DIR}/pass${PASS}/files_for_analysis/reads"
out_dir=${BASE_DIR}/pass${PASS}/files_for_analysis/fastqc/param_row_reads_${PARAM_ROW_READS}
mkdir -p $out_dir

# check if paired-end or single-end and call accordingly
if [[ "$paired_end_status" == "T" ]]; then
    echo "Detected paired-end reads"
    fqz_file1="${reads_dir}/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    fqz_file2="${reads_dir}/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz"
    fastqc -o $out_dir/ $fqz_file1 $fqz_file2
else
    echo "Detected single-end reads"
    fqz_file="${reads_dir}/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    fastqc -o $out_dir/ $fqz_file
fi


conda deactivate

module unload fastqc
