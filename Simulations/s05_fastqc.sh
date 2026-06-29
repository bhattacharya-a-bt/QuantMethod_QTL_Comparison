#!/bin/bash

module load fastqc

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate fastqc-0.11.9

PASS=$1
PARAM_ROW_READS=$2
GENO_PASS=$3
i=${LSB_JOBINDEX}
sample_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids"
sample_name=$(awk -v row=$i 'NR == row {print $1}' "$sample_file")

# path to your paramspace for reads file
psr_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/parameter_space_reads.txt"

# extract the paired_end status from the specified row and column
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

out_dir=/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/fastqc/param_row_reads_${PARAM_ROW_READS}
mkdir -p $out_dir

# check if paired-end or single-end and call accordingly
if [[ "$paired_end_status" == "T" ]]; then
    echo "Detected paired-end reads"
    fqz_file1="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    fqz_file2="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz"
    fastqc -o $out_dir/ $fqz_file1 $fqz_file2
else
    echo "Detected single-end reads"
    fqz_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
    fastqc -o $out_dir/ $fqz_file
fi


conda deactivate

module unload fastqc
