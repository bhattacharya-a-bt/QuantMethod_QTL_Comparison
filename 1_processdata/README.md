# Requantification of GTEx Files

For scripts related to the the processing and requantification of all GTEx v8 files for a given annotation and quantification method prior to downstream QTL/colocalization/TWAS, please see the `requantification_scripts` directory and its `README`.

# QTL / Colocalization / TWAS Pipeline

## Overview

Downstream of the requantification pipeline: takes the polished per-tissue, per-annotation, per-method expression RDS files and runs cis-eQTL mapping, GWAS colocalization, and TWAS model training/testing across all annotation × quantification-method × tissue combinations.

Most steps are driven by a shared 768 row parameter space file, `requant_paramspace.txt` (columns: `annot`, `quant`, `tissue`), with one LSF array index per row.

## s01 — Format covariates
**`r1_s01_writeCovariates.R`**
Loops over every tissue, reads GTEx's standard `.covar` covariate file, transposes/reformats it into the tab-delimited layout QTLtools expects, and writes `<tissue>_formatted_covariates.txt`. Run once per tissue, not part of the array jobs.

## s02 — RDS → BED
**`r1_s02_RDStoBed.R`** / **`r1_s02_RDStoBed.sh`**
For one annot/quant/tissue combo: loads the polished gene-level RDS, pulls out TMM-normalized counts, log2-transforms, filters low-expressed genes, attaches Ensembl gene coordinates, and writes a QTLtools-formatted BED of normalized expression (`<annot>_<quant>.v8.normalized_expression.bed`). Array job over the full paramspace (1-768).

## s03 — cis-eQTL mapping
**`r1_s03_qtltools.sh`** *(bash only, no R script)*
bgzip/tabix-indexes the BED from s02, then runs `QTLtools cis --permute 1000` against the GTEx genotype VCF and the formatted covariates to get permutation-pass cis-eQTL results per annot/quant/tissue. Array job (1–768).

*(s04 not included in this batch — presumably an eGene-aggregation step feeding the `r1_aggregated_eGene_lists.RDS` used in s05.)*

## s05 — Colocalization
**`r1_s05_coloc.R`** / **`r1_s05_coloc.sh`**
For a bin of eGenes (annot/quant/tissue + GWAS pheno file/name, with optional hg19→hg38 liftover): pulls GWAS summary stats near each gene, intersects SNPs with GTEx genotypes, computes LD, runs QTLtools nominal pass + eCAVIAR, and writes per-gene colocalization output (CLPP etc.) to `coloc_results/<tissue>/`.

## s06 — Train TWAS
**`r1_s06_gather_coloc.R`** / **`r1_s06_gather_coloc.sh`**
First, collects all per-gene coloc result files for one annot/quant/tissue across 7 phenotypes (breast/prostate cancer, T2D, BMI, height, schizophrenia, bipolar disorder), filters to significant colocalization hits, and writes one combined hits file per tissue/annot/quant/pheno. Array job (1–768).

**`r1_s06_trainTWAS.R`** / **`r1_s06_trainTWAS.sh`**
For one annot/quant/tissue + a gene-bin index: covariate-corrects the expression BED via `QTLtools correct`, then for each gene in the bin fits elastic net, BLUP, and SuSiE expression-prediction models (via `isotwas`) on local genotypes, keeps whichever has the best cross-validated R², and saves the SNP weights (or a placeholder if R² < 0.01) to `TWAS_weights/<tissue>_<annot>_<quant>/<gene>_TWAS.RDS`.

## s08 — Run TWAS
**`r1_s08_runTWAS.R`** / **`r1_s08_runTWAS.sh`**
For one row of the aggregated, R²-filtered TWAS model table (a tissue/annot/quant combo) and a bin of its genes: lifts over GWAS summary stats if needed, and for each gene intersects GWAS/model/genotype SNPs, computes LD, flips alleles as needed, and computes the TWAS burden Z-statistic plus the top local GWAS SNP. Writes one result file per gene/pheno/annot/quant.

## s09 — Gather TWAS results
**`r1_s09_gather_TWAS.sh`** 
For each annot/quant/tissue row in the paramspace and a fixed phenotype (`Height` as currently set), concatenates all per-gene TWAS Z-score files into one combined `r1_twas_z_<annot>_<quant>_<tissue>_<pheno>.txt`. Array job (1–768).

---

### Module/dependency notes
- R steps load a personal library path (`/rsrch5/home/epi/sthead/R/...`) and lean on `data.table`, `dplyr`, `edgeR`, `bigsnpr`, `rtracklayer`/`GenomicRanges` (liftover), and `isotwas` (TWAS model training).
- s03/s05/s06_trainTWAS/s08 all shell out to external binaries: `QTLtools`, `tabix`/`bgzip`, `plink2`, and (s05 only) `eCAVIAR`.
- Steps keyed on the paramspace file (`s02`, `s03`, `s06_gather_coloc`, `s09`) read their `annot`/`quant`/`tissue` triplet from `requant_paramspace.txt` via the LSF array index; `s05`, `s06_trainTWAS`, and `s08` take additional positional args (GWAS pheno info, gene bin index) on top of that.
