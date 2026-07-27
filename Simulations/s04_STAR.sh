#!/bin/sh

####################################################################################
# Usage:
#   bash run_star.sh <PASS> <PARAM_ROW_READS> <GENO_PASS> <THREADS> \
#                    <base_dir> <dir_temp> <txindex_star> <conda_bin>
#
# Arguments (positional, all required):
#   PASS            - pass number for reads / alignment / parameter files
#   PARAM_ROW_READS - row index into parameter_space_reads.txt
#   GENO_PASS       - pass number for the sample-IDs file
#   THREADS         - number of threads for STAR / samtools
#   base_dir        - base project directory (replaces the previously hardcoded root).
#                     Paths are reconstructed internally as:
#                       <base_dir>/pass<GENO_PASS>/files_for_analysis/1kg_eur_500_sample_ids
#                       <base_dir>/pass<PASS>/files_for_analysis/reads/
#                       <base_dir>/pass<PASS>/files_for_analysis/star_alignments/param_row_reads_<PARAM_ROW_READS>
#                       <base_dir>/pass<PASS>/files_for_analysis/parameter_space_reads.txt
#   dir_temp        - scratch directory for STAR temporary files
#                     (replaces the previously hardcoded bhattacharya_lab users path)
#   txindex_star    - full path to the STAR genome index directory
#                     (replaces the previously hardcoded GenomicReferences path)
#   conda_bin       - full path to the conda binary
#                     (e.g. /risapps/rhel8/miniforge3/24.5.0-0/bin/conda). Lives outside
#                     base_dir, so it is passed explicitly rather than reconstructed.
####################################################################################

# Set up paths and variables
PASS=$1
PARAM_ROW_READS=$2
GENO_PASS=$3
THREADS=$4
base_dir=$5
dir_temp=$6
txindex_star=$7
conda_bin=$8

i=${LSB_JOBINDEX}
sample_file="${base_dir}/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids"
SAMPLE=$(awk -v row=$i 'NR == row {print $1}' "$sample_file")

# directories for STAR alignment
DIR_TRIM="${base_dir}/pass${PASS}/files_for_analysis/reads" # Using FASTQ output as input for STAR
DIR_ALIGN="${base_dir}/pass${PASS}/files_for_analysis/star_alignments/param_row_reads_${PARAM_ROW_READS}"
BAM_FILE="${DIR_ALIGN}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

# path to paramspace for reads file
psr_file="${base_dir}/pass${PASS}/files_for_analysis/parameter_space_reads.txt"

# extract the paired_end status from the specified row and column
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

# Create output directories
mkdir -p ${DIR_ALIGN}
mkdir -p ${dir_temp}

eval "$(${conda_bin} shell.bash hook)"
conda activate samtools-1.16.1
module add star

# Run STAR alignment
echo "Running STAR alignment"
mkdir -p ${dir_temp}/${SAMPLE}
rm -rf ${dir_temp}/${SAMPLE}

# check if paired-end or single-end and call salmon accordingly
if [[ "$paired_end_status" == "T" ]]; then
  echo "Detected paired-end reads"
  STAR --genomeDir ${txindex_star} \
    --readFilesIn ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz \
    --runThreadN ${THREADS} \
    --readFilesCommand zcat \
    --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
    --sjdbScore 1 --outTmpDir ${dir_temp}/${SAMPLE} \
    --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 80000000000
else
  echo "Detected single-end reads"
  STAR --genomeDir ${txindex_star} \
    --readFilesIn ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz \
    --runThreadN ${THREADS} \
    --readFilesCommand zcat \
    --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
    --sjdbScore 1 --outTmpDir ${dir_temp}/${SAMPLE} \
    --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 80000000000
fi

samtools index ${BAM_FILE} -@ ${THREADS}

# Clean up temporary and intermediate files
rm -rf ${dir_temp}/${SAMPLE}
rm ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz
rm ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz

conda deactivate
