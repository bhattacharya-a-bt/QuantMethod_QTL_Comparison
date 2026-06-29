#!/bin/bash

tissue="$1"
user="$2"

cd /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts

bsub -env tissue=${tissue},user=${user} < autoLoop_process_short_reads_GR.lsf
