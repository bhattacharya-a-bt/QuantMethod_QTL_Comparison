#!/bin/bash
#BSUB -J "gatherColoc[1-768]"
#BSUB -o /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/coloc/gatherColoc_%J_%I.out
#BSUB -e /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/coloc/gatherColoc_%J_%I.err
#BSUB -q short
#BSUB -W 3:00
#BSUB -n 1
#BSUB -M 30
#BSUB -R rusage[mem=30]

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

# load modules

module load R/4.3.1

# assign analysis parameters

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

echo "Gathering colocalization results..."
Rscript /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts/r1_s06_gather_coloc.R ${annot} ${quant} ${tissue}

# unload modules
module unload R

echo "***** job ends ***** "
date

