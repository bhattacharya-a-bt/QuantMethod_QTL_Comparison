#!/bin/sh
# Master script for BAM to FASTQ conversion and Salmon quantification
# This script is called by the batch submission scripts

# Set up paths and variables
DIR_BAM=/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam
SAMPLES_FILE="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/bam_files_with_sample_attrib_rnaseq.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)
TISSUE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f9)
THREADS=5
OUT_DIR_SORTED="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/$TISSUE/sortedbam"
OUT_DIR_FASTQ="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/$TISSUE/fastq" # i also copied these over to /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/fastq

DIR_OUT_QUANT=/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/$TISSUE/salmon_quantifications/gencodev45
TXINDEX_salmon=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45

# Create output directories
mkdir -p ${OUT_DIR_SORTED}
mkdir -p ${OUT_DIR_FASTQ}
mkdir -p ${DIR_OUT_QUANT}

echo "Running task $LSB_JOBINDEX"
echo "LIBID: $LIBID"
echo "TISSUE: $TISSUE"
echo "SAMPLE: $SAMPLE"

# Activate conda environment
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1

# Step 1: Sort BAM file by read name
echo "Sorting BAM"
samtools sort -n -@ ${THREADS} -o ${OUT_DIR_SORTED}/${LIBID}_sorted.bam ${DIR_BAM}/${LIBID}

# Step 2: Convert BAM to FASTQ
echo "Dumping BAM to FASTQ"
samtools fastq -@ ${THREADS} -1 ${OUT_DIR_FASTQ}/${LIBID}_R1.fastq -2 ${OUT_DIR_FASTQ}/${LIBID}_R2.fastq -0 /dev/null -s /dev/null -n ${OUT_DIR_SORTED}/${LIBID}_sorted.bam

# Clean up sorted BAM file
rm ${OUT_DIR_SORTED}/${LIBID}_sorted.bam

# Step 3: Run Salmon quantification
salmon=/rsrch5/home/epi/bhattacharya_lab/software/salmon-latest_linux_x86_64/bin/salmon

echo "Running salmon"

$salmon quant \
  -i ${TXINDEX_salmon} --libType A \
  --validateMappings --seqBias \
  -p ${THREADS} \
  -1 ${OUT_DIR_FASTQ}/${LIBID}_R1.fastq \
  -2 ${OUT_DIR_FASTQ}/${LIBID}_R2.fastq \
  -o ${DIR_OUT_QUANT}/${SAMPLE}

# Clean up FASTQ files
rm ${OUT_DIR_FASTQ}/${LIBID}_R2.fastq
rm ${OUT_DIR_FASTQ}/${LIBID}_R1.fastq 

echo "Done with task $LSB_JOBINDEX"
echo "LIBID: $LIBID"
echo "TISSUE: $TISSUE"
echo "SAMPLE: $SAMPLE"

conda deactivate