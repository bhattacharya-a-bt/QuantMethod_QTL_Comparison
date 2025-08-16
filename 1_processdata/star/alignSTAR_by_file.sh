#!/bin/bash

# Usage: alignSTAR_by_file.sh <SAMPLES_FILE>
SAMPLES_FILE="$1"

if [ -z "$SAMPLES_FILE" ]; then
    echo "Usage: $0 <samples_file>"
    exit 1
fi

THREADS=16
DIR_BAM="/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam"
GENOME_DIR="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences"
TXINDEX_STAR="${GENOME_DIR}/star-2.7.4a_GCA_000001405.15_GRCh38_no_alt_analysis_set"

# Get total number of lines in file
TOTAL_LINES=$(wc -l < "$SAMPLES_FILE")
START_IDX=2  # skip header

for ((ROWID=START_IDX; ROWID<=TOTAL_LINES; ROWID++)); do
    LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
    SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)
    TISSUE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f9)

    echo "[$ROWID] Processing sample: $SAMPLE (Tissue: $TISSUE)"

    OUT_DIR_SORTED="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/sortedbam"
    OUT_DIR_FASTQ="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/fastq"
    DIR_ALIGN="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/star_alignments"
    DIR_TEMP="/rsrch5/scratch/epi/bhattacharya_lab/GTEx_reprocess/$TISSUE/star_temp"

    BAM_FILE="${DIR_ALIGN}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    BAM_INDEX_FILE="${BAM_FILE}.bai"
    SJ_TAB_FILE="${DIR_ALIGN}/${SAMPLE}_SJ.out.tab"

    if [ ! -s "$BAM_FILE" ] || [ ! -s "$BAM_INDEX_FILE" ] || [ ! -s "$SJ_TAB_FILE" ] || [ "$(stat -c%s "$BAM_FILE")" -le 1000000000 ]; then
        echo "Files missing or too small. Proceeding with alignment for $SAMPLE."

        eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
        conda activate samtools-1.16.1
        module add star

        mkdir -p "$OUT_DIR_SORTED" "$OUT_DIR_FASTQ" "$DIR_ALIGN" "$DIR_TEMP"

        echo "Sorting BAM for $SAMPLE"
        samtools sort -n -@ $THREADS -o ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam ${DIR_BAM}/${LIBID}

        echo "Dumping FASTQ for $SAMPLE"
        samtools fastq -@ $THREADS -1 ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq -2 ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq -0 /dev/null -s /dev/null -n ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam
        rm -f ${OUT_DIR_SORTED}/${SAMPLE}_sorted.bam

        echo "Running STAR alignment for $SAMPLE"
        rm -rf ${DIR_TEMP}/${SAMPLE}

        STAR --genomeDir ${TXINDEX_STAR} \
            --readFilesIn ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq \
            --runThreadN ${THREADS} \
            --genomeLoad NoSharedMemory \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMheaderHD @HD VN:1.4 SO:coordinate \
            --outSAMunmapped Within \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbScore 1 \
            --outTmpDir ${DIR_TEMP}/${SAMPLE} \
            --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
            --outBAMsortingBinsN 200 \
            --limitBAMsortRAM 80000000000

        samtools index ${BAM_FILE} -@ ${THREADS}

        echo "Cleaning up temp and FASTQ for $SAMPLE"
        rm -rf ${DIR_TEMP}/${SAMPLE}
        rm -f ${OUT_DIR_FASTQ}/${SAMPLE}_R1.fastq ${OUT_DIR_FASTQ}/${SAMPLE}_R2.fastq

        conda deactivate
    else
        echo "All STAR output files exist and are of sufficient size. Skipping $SAMPLE."
    fi
done
