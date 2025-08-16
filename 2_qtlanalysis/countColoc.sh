#!/bin/bash
#BSUB -J countcoloc[1-480]
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/countcoloc_%J_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/countcoloc_%J_%I.err
#BSUB -W 24:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q medium

# Load required modules (if needed)
module load R/4.4.3

# Calculate tissue_index and cancer_index from LSB_JOBINDEX
# We have 48 tissues and 10 cancers = 480 total combinations
# Job indices 1-480 map to:
# Job 1: tissue=1, cancer=1
# Job 2: tissue=1, cancer=2
# ...
# Job 10: tissue=1, cancer=10
# Job 11: tissue=2, cancer=1
# ...
# Job 480: tissue=48, cancer=10

tissue_index=$(((LSB_JOBINDEX - 1) / 10 + 1))
cancer_index=$(((LSB_JOBINDEX - 1) % 10 + 1))

echo "Job $LSB_JOBINDEX: Processing tissue_index=$tissue_index, cancer_index=$cancer_index"

# Run the R script with the calculated indices
cd /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis
Rscript countColoc.R $tissue_index $cancer_index