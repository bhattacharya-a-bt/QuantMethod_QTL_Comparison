#!/bin/sh

# Set up paths and variables
PASS=$1
PARAM_ROW_READS=$2
GENO_PASS=$3
THREADS=$4

i=${LSB_JOBINDEX}
sample_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids"
SAMPLE=$(awk -v row=$i 'NR == row {print $1}' "$sample_file")

# directories for STAR alignment
DIR_TRIM="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/reads"  # Using FASTQ output as input for STAR
DIR_ALIGN="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/star_alignments/param_row_reads_${PARAM_ROW_READS}"
DIR_TEMP="/rsrch5/scratch/epi/bhattacharya_lab/users/sthead"
BAM_FILE="${DIR_ALIGN}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

# STAR genome reference
GENOME_DIR="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences"
TXINDEX_STAR="${GENOME_DIR}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set"

# path to paramspace for reads file
psr_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/parameter_space_reads.txt"

# extract the paired_end status from the specified row and column
paired_end_status=$(awk -v row=$((PARAM_ROW_READS + 1)) 'NR == row {print $3}' "$psr_file")

# Create output directories
mkdir -p ${DIR_ALIGN}
mkdir -p ${DIR_TEMP}

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1
module add star

# Run STAR alignment
echo "Running STAR alignment"
mkdir -p ${DIR_TEMP}/${SAMPLE}
rm -rf ${DIR_TEMP}/${SAMPLE}

# check if paired-end or single-end and call salmon accordingly
if [[ "$paired_end_status" == "T" ]]; then
    echo "Detected paired-end reads"
    STAR --genomeDir ${TXINDEX_STAR} \
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
    --sjdbScore 1 --outTmpDir ${DIR_TEMP}/${SAMPLE} \
    --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 80000000000
else
    echo "Detected single-end reads"
    STAR --genomeDir ${TXINDEX_STAR} \
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
    --sjdbScore 1 --outTmpDir ${DIR_TEMP}/${SAMPLE} \
    --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 80000000000
fi




samtools index ${BAM_FILE} -@ ${THREADS}

# Clean up temporary and intermediate files
rm -rf ${DIR_TEMP}/${SAMPLE}
rm ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R1.fastq.gz 
rm ${DIR_TRIM}/sim_${SAMPLE}_param_row_reads_${PARAM_ROW_READS}_R2.fastq.gz

conda deactivate

