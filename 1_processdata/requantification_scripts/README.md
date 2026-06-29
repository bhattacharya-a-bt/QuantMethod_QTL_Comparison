# GTEx v8 Requantification Scripts

Scripts for regenerating GTEx v8 RNA-seq quantifications from the original BAM files, across **4 transcript annotation versions** and **4 quantification methods**, for the quantification-method/annotation TWAS comparison project. Runs on an LSF cluster.

**Annotations:** GENCODE v27 (GTEx v8's default), GENCODE v38, GENCODE v45, Ensembl v115 (all GRCh38)
**Methods:** Salmon (alignment-free), STAR+featureCounts, STAR+RSEM, kallisto
**Levels:** gene and transcript, with raw counts + TMM-normalized CPM for each

## Pipeline at a glance

```
1. Inventory source BAMs           count_gtex_bams.sh, check_GTEx_files.sh, check_GTEx_files.R
2. Per-sample align + quantify      loop_process_short_reads_GR.sh
                                       -> autoLoop_process_short_reads_GR.lsf  (disk-throttled job submitter)
                                            -> run_process_short_reads_GR.lsf  (the actual per-sample worker)
                                     loop10_process_short_reads_GR.lsf         (10-sample test harness)
3. Aggregate to per-tissue RDS      run_quants_GR.sh -> run_quants_GR.lsf -> quants_GR.R
4. QC + harmonize across methods    check_outputs.R, polish_outputs.R, count_polished_samples.R
5. Free disk                        cleanup_GR.sh, cleanup_GR_rsrch8.sh
```

Steps 2 and 3 are independent passes (per-sample quant, then per-tissue aggregation) and are typically run tissue-by-tissue: submit step 2 for a tissue, wait for all samples to finish, then submit step 3 for that tissue.

---

## 1. Inventory scripts

**`count_gtex_bams.sh <tissue>`**
Counts how many sample IDs for a given tissue (from the GTEx sample attributes file) actually have a matching BAM in the source BAM directory. Quick sanity check before launching a full run.

**`check_GTEx_files.sh`**
Walks every sample ID in the GTEx sample attributes file and writes a `sample_id -> bam_file` lookup table (`GTE_bam_list.tsv`) by globbing the source BAM directory.

**`check_GTEx_files.R`**
Cross-checks `GTE_bam_list.tsv` against GTEx's own `sequencing.json` manifest (filtered to `Aligned Reads` / `bam` / `Sequencing Reads`), joins in tissue labels from the sample attributes file, and prints a per-tissue summary of how many BAMs were found vs. expected.

## 2. Per-sample alignment and quantification

**`run_process_short_reads_GR.lsf`** — the worker job (12 cores, 64 GB, `medium` queue, 24 h limit). Takes a single sample via environment variables (`FILE`, `TISSUE`, `USER`, plus boolean flags `PULL`/`ALIGN`/`SALMON`/`FEATURECOUNTS`/`RSEM`/`KALLISTO`) and, depending on which flags are set:
- `PULL`: extracts paired FASTQs from the source BAM with `samtools sort -n` + `samtools fastq`
- `ALIGN`: 2-pass STAR alignment to GRCh38 (needed only for featureCounts)
- `FEATURECOUNTS` / `RSEM` / `KALLISTO` / `SALMON`: each, if set, runs that method across **all four annotation versions in a single loop** (separate index/annotation per version), then deletes its temporary inputs (STAR BAM for featureCounts, RSEM's intermediate BAM, etc.)
- Cleans up the pulled FASTQs at the end if `PULL` was used

Conda environments used: `samtools-1.16.1`, `star-2.7.4a`, `rsem-1.3.3`, `kallisto-0.48.0`, `salmon-1.10.2`, plus an environment module `subread` for featureCounts.

**`autoLoop_process_short_reads_GR.lsf`** — orchestrator (1 core, `vlong` queue, 7-day limit) that drives `run_process_short_reads_GR.lsf` across an entire tissue:
1. Lists all sample IDs for the tissue.
2. For each sample, checks whether its **Ensembl** outputs exist for each method (`_transcript.txt`, `.isoforms.results`, `quant.sf`, `abundance.tsv`) to decide which `SALMON`/`FEATURECOUNTS`/`RSEM`/`KALLISTO` flags are still needed — Ensembl is checked as a stand-in for "all four annotations," since one method-run inside the worker job covers all four versions together. Samples where everything already exists are skipped.
3. Throttles concurrent job submission by **free scratch disk space** rather than a fixed job count: estimates 80 GB per active sample plus a 500 GB safety buffer, checks `bjobs` every 5 minutes, and submits more samples only as disk frees up.
4. Waits for all jobs to finish before exiting.

**`loop_process_short_reads_GR.sh <tissue> <user>`** — the actual entry point: just `bsub`s `autoLoop_process_short_reads_GR.lsf` with the right environment variables.

**`loop10_process_short_reads_GR.lsf`** — a debug/validation variant of the autoloop with the tissue (`Brain_Hippocampus`) and user hard-coded. Selects only the first 10 samples that have a BAM, submits them without disk throttling, then runs the same 5-minute disk-usage monitor loop. Intended for a smoke test before committing a full tissue to the queue — **edit the hard-coded `tissue`/`user` at the top before reusing.**

## 3. Aggregation to per-tissue R objects

**`run_quants_GR.sh <tissue> <user>`** — wrapper that `bsub`s `run_quants_GR.lsf`.

**`run_quants_GR.lsf`** — LSF job (12 cores, 200 GB, `medium` queue, 8 h limit) that runs `quants_GR.R` inside a Singularity container (`rstudio_4.3.1.sif`).

**`quants_GR.R`** — the core aggregation script. For one tissue, loops over all four annotation versions and, within each, runs four independent (each wrapped in `tryCatch`, so one failure doesn't kill the run) import-and-save steps:

| Method | Import path | Gene-level object | Transcript-level object |
|---|---|---|---|
| Salmon | `tximeta` on each sample's `quant.sf` | `SummarizedExperiment` (via `summarizeToGene`) | `SummarizedExperiment` |
| STAR + featureCounts | parallel `fread` + merge of per-sample `_gene.txt`/`_transcript.txt` count files (SAF-based) | `list(counts, TMM_normalized)` | `list(counts, TMM_normalized)` |
| STAR + RSEM | `tximport(type="rsem")` on `.genes.results`/`.isoforms.results` | `list(counts, TMM_normalized)` | `list(counts, TMM_normalized)` |
| kallisto | `tximport(type="kallisto")` on `abundance.tsv`, summarized via tx2gene map | `SummarizedExperiment` | `SummarizedExperiment` |

Other notable behavior:
- Each method/level combination is **skipped if its output RDS already exists**, so the script is safely resumable / re-runnable.
- Before importing, each method does its own QC pass: drops samples with missing/empty/truncated quant files (and, for RSEM, drops any sample whose gene/transcript count doesn't match the most common count across samples — a proxy for a corrupted or mismatched annotation run).
- TMM normalization (`edgeR::calcNormFactors` + `cpm`) is applied at both gene and transcript level. For >100 samples, normalization factors are computed once on the **full** sample set and then CPM is computed in chunks of 100 samples at a time (to bound memory) — chunking never changes the factors used, only how the matrix multiplication is batched.
- All errors and warnings across the whole run are collected and written to `<tissue>_combine.log` in the requants directory, alongside a console summary of error/warning counts.
- Output directory: `requants/<ANNOTATION>/`, e.g. `requants/GENCODE_v38/Whole_Blood_salmon_gene.RDS`.

## 4. QC and cross-method/annotation harmonization

**`check_outputs.R`** — scans all four annotation directories for `<tissue>_<method>_<level>.RDS` files, infers the complete tissue list from whatever's present, and reports (to console and to `file_check_report.txt`) exactly which of the 4 annotations × 4 methods × 2 levels = 32 expected files per tissue are missing, with an overall completion percentage.

**`polish_outputs.R`** — because each method/annotation combination can independently drop different samples or features during QC, this script finds the common ground so that downstream eQTL/TWAS comparisons are apples-to-apples:
1. For each tissue, for each annotation, intersects sample IDs across all 4 methods, and intersects gene IDs (and transcript IDs) across all 4 methods.
2. Intersects sample IDs again across all 4 annotations, so the **final sample set is identical for every method and every annotation** for that tissue.
3. Subsets every RDS file to the common samples/features and writes the result under `requants/polished/<ANNOTATION>/`.
4. Skips a tissue entirely if any of its 32 input files are missing, and skips re-polishing if the polished outputs already exist.

**`count_polished_samples.R`** — quick sanity check that loads each tissue's polished `*_featureCounts_gene.RDS` (Ensembl only, as currently written) and reports the sample count, to spot tissues with unexpectedly small final sample sizes after harmonization.

## 5. Cleanup

**`cleanup_GR.sh <tissue> <user>`** — once a tissue's RDS aggregates exist, this deletes the bulky per-sample intermediates: scratch-space FASTQs/alignments/temp dirs, and the raw per-sample quant subdirectories (`salmon/`, `featureCounts/`, `RSEM/`, `kallisto/`) under each of the four annotation directories. 

## Dependencies

- **Scheduler:** LSF (`bsub`, `bjobs`)
- **Conda envs:** `samtools-1.16.1`, `star-2.7.4a`, `rsem-1.3.3`, `kallisto-0.48.0`, `salmon-1.10.2`
- **R 4.3.1** (via Singularity image `rstudio_4.3.1.sif`), packages: `tximeta`, `tximport`, `data.table`, `stringr`, `edgeR`, `SummarizedExperiment`, `Rsubread`, `rtracklayer`, `pbapply`, `parallel`, `jsonlite`, `dplyr`, `tidyr`
- Reference genome: GRCh38 (`GCA_000001405.15_no_alt_analysis_set`); per-annotation Salmon/kallisto indices, RSEM references, and featureCounts SAF files are assumed to already exist on disk.


