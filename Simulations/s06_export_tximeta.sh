#!/bin/bash

# load modules
module load R/4.3.1

PASS=$1
PARAM_ROW_READS=$2
ANNOT=$3

# run the simulationR script with the job index as an argument
Rscript /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/scripts/s06_export_tximeta.R ${PASS} ${PARAM_ROW_READS} ${ANNOT}

# unload modules
module unload R

echo "***** job ends ***** "
date
