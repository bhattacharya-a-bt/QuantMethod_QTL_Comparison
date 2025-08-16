#!/bin/bash
#BSUB -J gtextwas[1-576]
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/twas_%J_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/twas_%J_%I.err
#BSUB -W 24:00
#BSUB -n 10
#BSUB -M 100
#BSUB -R rusage[mem=80]
#BSUB -q medium


# Load modules or activate R environment if needed
module load plink htslib tabix qtltools
module load R/4.4.3

# List of tissues (modify this if not alphabetically ordered)
TISSUE_LIST=( $(ls /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/gencodev45) )

# Calculate tissue index and secondary index
# Each tissue gets 12 jobs (indices 1-12)
TISSUE_INDEX=$(( (LSB_JOBINDEX - 1) / 12 ))
SECONDARY_INDEX=$(( ((LSB_JOBINDEX - 1) % 12) + 1 ))

TISSUE=${TISSUE_LIST[$TISSUE_INDEX]}

echo "Running job for tissue: $TISSUE"
echo "Secondary index: $SECONDARY_INDEX"
echo "Job array index: $LSB_JOBINDEX"

cd /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas
Rscript trainTWAS.R $TISSUE $SECONDARY_INDEX