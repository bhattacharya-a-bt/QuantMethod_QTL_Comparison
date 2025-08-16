#!/bin/bash
#BSUB -J star_align[48,96]
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/1_processdata/star/logs/star_align_%J_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/1_processdata/star/logs/star_align_%J_%I.err
#BSUB -W 200:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q long

# Set environment
TISSUE_DIR="/rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/tissue_samples"
SCRIPT="/rsrch5/home/epi/abhattacharya3/processGTEx/1_processdata/star/alignSTAR_by_file.sh"

# Make sure tissue list is in consistent order
TISSUE_LIST=$(ls "$TISSUE_DIR"/*.txt | sort)

# Calculate index and version based on LSB_JOBINDEX
# Jobs 1-48: index 1-48, version 38
# Jobs 49-96: index 1-48, version 45
if [ $LSB_JOBINDEX -le 48 ]; then
    INDEX=$LSB_JOBINDEX
    VERSION=38
else
    INDEX=$((LSB_JOBINDEX - 48))
    VERSION=45
fi

TISSUE_FILE=$(echo "$TISSUE_LIST" | sed -n "${INDEX}p")

echo "Running job $LSB_JOBINDEX: index=$INDEX, version=$VERSION"
echo "Processing tissue file: $TISSUE_FILE"

# Load R module and run script
module load R/4.4.3
Rscript /rsrch5/home/epi/abhattacharya3/processGTEx/1_processdata/star/countGene_featureCounts.R --index $INDEX --version $VERSION