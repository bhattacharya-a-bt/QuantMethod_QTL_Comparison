#!/bin/sh
# Master script for BAM to FASTQ conversion and STAR alignment
# This script is called by the batch submission scripts

# Set up paths and variables
DIR_BAM="/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam"
SAMPLES_FILE="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/bam_files_with_sample_attrib_rnaseq.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)
TISSUE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f9)
THREADS=8

echo "Running task $LSB_JOBINDEX"
echo "LIBID: $LIBID"
echo "TISSUE: $TISSUE"
echo "SAMPLE: $SAMPLE"

# Original output directories
OUT_DIR_SORTED="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/sortedbam"
OUT_DIR_FASTQ="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/fastq"

# New directories for STAR alignment
DIR_TRIM="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/fastq"  # Using FASTQ output as input for STAR
DIR_ALIGN="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/star_alignments"
DIR_TEMP="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/star_temp"

# STAR genome reference
GENOME_DIR=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences
TXINDEX_STAR=${GENOME_DIR}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set

# Create output directories
mkdir -p ${OUT_DIR_SORTED}
mkdir -p ${OUT_DIR_FASTQ}
mkdir -p ${DIR_ALIGN}
mkdir -p ${DIR_TEMP}

# Output file paths
BAM_FILE="${DIR_ALIGN}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
BAM_INDEX_FILE="${BAM_FILE}.bai"
SJ_TAB_FILE="${DIR_ALIGN}/${SAMPLE}_SJ.out.tab"

# Check if all expected files exist and are valid
if [ ! -s "$BAM_FILE" ] || [ ! -s "$BAM_INDEX_FILE" ] || [ ! -s "$SJ_TAB_FILE" ] || [ "$(stat -c%s "$BAM_FILE")" -le 1000000000 ]; then
    echo "Files missing or too small. Proceeding with alignment for $SAMPLE."

    # Activate conda environment
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate samtools-1.16.1
    module add star

    # Step 1: Sort BAM file by read name
    echo "Sorting BAM"
    rm -f ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam
    samtools sort -n -@ ${THREADS} -o ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam ${DIR_BAM}/${LIBID}

    # Step 2: Convert BAM to FASTQ
    echo "Dumping BAM to FASTQ"
    rm -f ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq
    samtools fastq -@ ${THREADS} -1 ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq -2 ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq -0 /dev/null -s /dev/null -n ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam

    # Clean up sorted BAM file
    rm -f ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam

    # Step 3: Run STAR alignment
    echo "Running STAR alignment"
    mkdir -p ${DIR_TEMP}/${SAMPLE}
    rm -rf ${DIR_TEMP}/${SAMPLE}

    STAR --genomeDir ${TXINDEX_STAR} \
        --readFilesIn ${DIR_TRIM}/${SAMPLE}_R1.fastq ${DIR_TRIM}/${SAMPLE}_R2.fastq \
        --runThreadN ${THREADS} \
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

    samtools index ${BAM_FILE} -@ ${THREADS}

    # Clean up temporary and intermediate files
    rm -rf ${DIR_TEMP}/${SAMPLE}
    rm -f ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq

    conda deactivate

else

    echo "All STAR output files exist and are of sufficient size. Skipping alignment for $SAMPLE."
fi
