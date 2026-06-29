#!/bin/bash
#BSUB -J "make_bed[273-288]"
#BSUB -o /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/qtl/make_bed_%J_%I.out
#BSUB -e /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/qtl/make_bed_%J_%I.out
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 20
#BSUB -R rusage[mem=20]

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "User: $LSB_USER"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

# load modules
module load R/4.3.1

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

Rscript /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts/r1_s02_RDStoBed.R ${annot} ${quant} ${tissue}

# unload modules
module unload R

echo "***** job ends ***** "
date
