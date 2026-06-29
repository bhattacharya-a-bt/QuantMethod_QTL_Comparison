#!/bin/bash

# load modules
module load R/4.3.1

PASS=$1
PARAM_ROW_READS=$2

# run the simulationR script with the job index as an argument
Rscript /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/scripts/s02_sim_reads.R ${LSB_JOBINDEX} ${PASS} ${PARAM_ROW_READS}

# unload modules
module unload R

echo "***** job ends ***** "
date
