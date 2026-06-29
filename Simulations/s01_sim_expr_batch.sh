#!/bin/bash

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
module load tabix/0.2.6

PASS=$1
GENO_PASS=$2
GENES_PER_JOB=$3
TOTAL_GENES=$4
TOTAL_JOBS=$5

JOB_INDEX=${LSB_JOBINDEX}
START_INDEX=$(( (JOB_INDEX - 1) * GENES_PER_JOB + 1 ))
END_INDEX=$(( JOB_INDEX * GENES_PER_JOB ))

if [ $END_INDEX -gt $TOTAL_GENES ]; then
    END_INDEX=$TOTAL_GENES
fi

for ((i=START_INDEX; i<=END_INDEX; i++))
do
    Rscript /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/scripts/s01_sim_expr.R $i ${PASS} ${GENO_PASS}
done

# unload modules
module unload R
module unload tabix

echo "***** job ends ***** "
date