#!/bin/bash
#BSUB -J "clean_qtl_[1-22]"
#BSUB -o /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/logs/qtl/get_true_beta_info_%I.out
#BSUB -e /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/logs/qtl/get_true_beta_info_%I.err
#BSUB -n 1
#BSUB -M 30
#BSUB -R "rusage[mem=30G]"
#BSUB -W 24:00
#BSUB -q medium

# load modules
module load R/4.3.1

CHR=$LSB_JOBINDEX

# run the simulationR script with the job index as an argument
Rscript /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/scripts/s10_clean_qtl_res.R ${CHR}

# unload modules
module unload R

echo "***** job ends ***** "
date
