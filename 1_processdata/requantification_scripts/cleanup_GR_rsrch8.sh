#!/bin/bash

tissue="$1"
user="$2"

rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/temp/${tissue}
rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/raw/${tissue}
rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/reports/${tissue}
rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/metrics/${tissue}
rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/align/${tissue}
rm -rf /rsrch8/scratch/epi/${user}/GTEx_v8/quant/${tissue}
rm -rf /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/GENCODE_v27/${tissue}
rm -rf /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/GENCODE_v38/${tissue}
rm -rf /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/GENCODE_v45/${tissue}
rm -rf /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Ensembl/${tissue}