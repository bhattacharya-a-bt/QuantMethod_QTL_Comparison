#!/bin/bash

tissue="$1"
user="$2"

DIR_SCRIPTS=/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts

bsub -env TISSUE=${tissue},USER=${user} < ${DIR_SCRIPTS}/run_quants_GR.lsf
