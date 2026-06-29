# Simulation Pipeline

## Quick Overview

Simulates RNA-seq data with a known ground-truth isoform-level eQTL architecture, then runs it through a subset of the quant/eQTL machinery as the GTEx pipeline. Organized into numbered "passes" (e.g. `pass1`, `pass2`) under `GTEx_gencode_comp/pass<N>/`; `geno_pass` lets a later pass reuse an earlier pass's genotype subset instead of regenerating it.

## Setup (run once per pass / genotype pass)

**`prior_s00_set_parameter_space_pass2.R`**
Picks the genes/isoforms to simulate: filters GENCODE v38 to autosomal protein-coding genes with ≥1 isoform, optionally can subsample to a fixed number of genes, and for each gene draws a true heritability (h²) and an isoform-sharing parameter (how many causal eQTLs are shared vs. isoform-specific). Writes `gene_level_parameters.txt` (the per-gene "answer key") and `parameter_space_reads.txt` (library size / read length / paired-end combos to simulate).

**`prior_s01_subset_1KG.R`**
Randomly selects 500 European 1000 Genomes samples to use as the simulated genotype cohort; writes the sample ID list.

**`prior_s02_convert_and_filt_vcf.sh`**
For each chromosome, subsets the 1KG hg38 PLINK files to those 500 samples, filters to MAF ≥ 0.05 SNPs, and exports/bgzips/tabixes a per-chromosome VCF.

## s01 — Simulate isoform expression
**`s01_sim_expr.R`** / **`s01_sim_expr_batch.sh`**
For one gene: pulls local genotypes from the 1KG VCF, simulates a multivariate isoform-expression model with shared + isoform-specific causal eQTLs (effect sizes drawn to hit the gene's target h² and an isoform-correlation structure), and saves the simulated expression matrix plus the true SNP effects/heritability (`B`, `isoform_h2`, `gene_h2`). The batch script chunks genes across an LSF array.

## s02 — Simulate reads
**`s02_sim_reads.R`** / **`s02_sim_reads.sh`**
Aggregates all per-gene simulated expression into one matrix, converts one sample's simulated (log) expression to linear scale, and calls `Rsubread::simReads` against the GENCODE v38 transcriptome FASTA to generate that sample's simulated FASTQ(s) at the specified library size/read length/paired-end setting.

## s03–s04 — Quantify / align simulated reads
**`s03_salmon.sh`** — runs Salmon (`quant`, `--validateMappings --seqBias`) on the simulated reads against a given transcriptome index, per sample.
**`s04_STAR.sh`** — aligns the simulated reads to the genome with STAR (coordinate-sorted BAM, standard ENCODE-style params), indexes with samtools, and cleans up temp files/input FASTQs.

## s05 — QC + gene-level counts
**`s05_fastqc.sh`** → **`s06_multiqc.sh`** + **`s06_multiqc.R`** — FastQC per sample, MultiQC to roll the per-sample reports into one summary table, then the R script flags any sample that failed *any* FastQC module and writes the passing-sample list.
**`s05_featureCounts.R`** / **`s05_featureCounts.sh`** — runs `Rsubread::featureCounts` on the STAR BAMs against a given GTF/annotation, saving raw gene-level counts for that annotation.

## s06 — Aggregate quantifications
**`s06_export_tximeta.R`** / **`s06_export_tximeta.sh`**
Imports the Salmon quant.sf files via `tximeta` (transcript + gene level) and the featureCounts RDS from s05, TMM-normalizes both (plus a log2(TMM+1) assay), and saves everything as `SummarizedExperiment` RDS objects.

## s08 — Build QTLtools input BEDs
**`s08_prep_qtl_files.R`**
For each of 4 annotation versions (`gencode_v27/v38/v45`) and both methods (Salmon, featureCounts): computes TPM, performs expression QC, and writes QTLtools-formatted BED files for each of three considered value types (raw abundance, TMM-normalized, log2-TMM) — one BED per annotation × method × value-type combination.

## s09 — cis-eQTL mapping
**`s09_qtltools.sh`**
Per chromosome: bgzip/tabix-indexes the matching BED from s08, then runs `QTLtools cis --permute 1000` against the per-chromosome simulated genotype VCF from `prior_s02`, for a given annotation/measure/method, optionally rank-normalizing (`--normal`).

## s10 — Validate against ground truth
**`s10_clean_qtl_res.R`** / **`s10_clean_qtl_res.sh`**
Per chromosome: loads the QTLtools eQTL results plus the true per-gene simulation objects from `s01` (`B`, `isoform_h2`), matches detected variants to the true causal SNPs (via rsID lookup against the 1KG `.bim` and LD with `plink2 --r2-unphased`), and flags each called eQTL as a true/shared/isoform-specific hit or not.

---

### Notes
- "Pass" numbers track simulation iterations (parameter tweaks, gene-set changes, etc.); `geno_pass` lets you point at an earlier pass's 1KG genotype subset/VCFs without redoing `prior_s01`/`prior_s02`.
- This simulation arm only compares **Salmon vs. featureCounts** (not RSEM/kallisto, which appear in the applied GTEx requantification scripts) across the same three GENCODE versions. We further do not consider Ensembl here in simulations.
- `param_row_reads` indexes into `parameter_space_reads.txt` (library size × read length × paired-end) — there's currently just one row defined, but the scripts are written to support a full grid.

