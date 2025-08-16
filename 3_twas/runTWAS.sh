#!/bin/bash
#BSUB -J "twas_analysis[1-1000]"
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/twas_analysis_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs/twas_analysis_%I.err
#BSUB -n 1
#BSUB -M 40
#BSUB -R rusage[mem=40]
#BSUB -W 96:00
#BSUB -q long

# Load required modules (adjust as needed for your system)
module load R/4.4.3 plink

# Calculate array_index and gwas_index from LSB_JOBINDEX
# We have 100 array chunks and 10 GWAS files = 1000 total jobs
# Job 1-100: array_index=1-100, gwas_index=1
# Job 101-200: array_index=1-100, gwas_index=2
# etc.

gwas_index=$(( (${LSB_JOBINDEX} - 1) / 100 + 1 ))
array_index=$(( (${LSB_JOBINDEX} - 1) % 100 + 1 ))

echo "Job ${LSB_JOBINDEX}: Processing array_index=${array_index}, gwas_index=${gwas_index}"

# Create logs directory if it doesn't exist
mkdir -p /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/logs

# Run the R script with the calculated indices
Rscript /rsrch5/home/epi/abhattacharya3/processGTEx/4_twas/runTWAS.R ${array_index} ${gwas_index}