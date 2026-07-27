#!/bin/bash

####################################################################################
# Usage:
#   bash run_salmon.sh <PASS> <PARAM_ROW_READS> <GENO_PASS> <TXINDEX_salmon> \
#                      <TXOME_name> <THREADS> <base_dir> <salmon>
#
# Arguments (positional, all required):
#   PASS            - pass number for reads / parameter / output files
#   PARAM_ROW_READS - row index into parameter_space_reads.txt
#   GENO_PASS       - pass number for the sample-IDs file
#   TXINDEX_salmon  - path to the salmon transcriptome index
#   TXOME_name      - transcriptome name (used in output dir label)
#   THREADS         - number of threads for salmon
#   base_dir        - base project directory (replaces the previously hardcoded root).
#                     Paths are reconstructed internally as:
#                       <base_dir>/pass<PASS>/files_for_analysis/        (working dir)
#                       <base_dir>/pass<PASS>/files_for_analysis/reads/
#                       <base_dir>/pass<PASS>/files_for_analysis/parameter_space_reads.txt
#                       <base_dir>/pass<GENO_PASS>/files_for_analysis/1kg_eur_500_sample_ids
#   salmon          - full path to the salmon binary
#                     (e.g. .../salmon-latest_linux_x86_64/bin/salmon). Lives outside
#                     base_dir, so it is passed explicitly rather than reconstructed.
####################################################################################

PASS=$1
PARAM_ROW_READS=$2
GENO_PASS=$3
TXINDEX_salmon=$4
TXOME_name=$5
THREADS=$6
base_dir=$7
salmon=$8

cd ${base_dir}/pass${PASS}/files_for_analysis

filename="reads_file_list_${PARAM_ROW_READS}"
# check if the file does not exist
if [ ! -f "$filename" ]; then
  # create the file
  ls reads/sim*_${PARAM_ROW_READS}_* > $filename
  echo "$filename created."
else
  echo "$filename already exists."
fi

# path to paramspace for reads file
psr_file="${base_dir}/pass${PASS}/files_for_analysis/parameter_space_reads.txt"

# extract the paired_end status from the specified row and column
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

i=${LSB_JOBINDEX}
sample_file="${base_dir}/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids"
sample_name=$(awk -v row=$i 'NR == row {print $1}' "$sample_file")

out_dir="salmon/param_row_reads_${PARAM_ROW_READS}_${TXOME_name}"

mkdir -p "$out_dir"

# have no dumping of EC data at the moment to save memory

# check if paired-end or single-end and call salmon accordingly
if [[ "$paired_end_status" == "T" ]]; then
  echo "Detected paired-end reads"
  fqz_file1="${base_dir}/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
  fqz_file2="${base_dir}/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz"
  # Call salmon for paired-end reads
  $salmon quant -i ${TXINDEX_salmon} -l A \
    -1 ${fqz_file1} -2 ${fqz_file2} \
    -o ${out_dir}/${sample_name} \
    -p ${THREADS} \
    --validateMappings \
    --seqBias
else
  echo "Detected single-end reads"
  fqz_file="${base_dir}/pass${PASS}/files_for_analysis/reads/sim_${sample_name}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz"
  # Call salmon for single-end reads
  $salmon quant -i ${TXINDEX_salmon} -l A \
    --unmatedReads ${fqz_file} \
    -o ${out_dir}/${sample_name} \
    -p ${THREADS} \
    --validateMappings \
    --seqBias
fi

module unload salmon
