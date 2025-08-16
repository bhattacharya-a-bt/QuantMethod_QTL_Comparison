#!/bin/bash
#BSUB -J "coloc_analysis[1-1536]"
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/coloc_analysis_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/coloc_analysis_%I.err
#BSUB -n 1
#BSUB -M 40
#BSUB -R rusage[mem=40]
#BSUB -W 48:00
#BSUB -q long

# Load required modules (adjust as needed for your system)
module load R/4.4.3 plink

# Run the R script with the array index
Rscript /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/runColoc.R ${LSB_JOBINDEX}