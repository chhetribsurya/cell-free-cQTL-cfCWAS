# Roadmap Bootstrap and Leave-One-Out Overlap

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![BEDTools](https://img.shields.io/badge/BEDTools-2.30%2B-green.svg)](https://bedtools.readthedocs.io/)
[![Bash](https://img.shields.io/badge/Bash-4.0%2B-lightgrey.svg)](https://www.gnu.org/software/bash/)
[![Roadmap](https://img.shields.io/badge/Roadmap-18--state%20ChromHMM-orange.svg)](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html)

This directory contains two robustness analyses for cfChIP overlap with **Roadmap / EpiMap 18-state ChromHMM** reference epigenomes (H3K27ac-dense tracks):

1. **Bootstrap overlap** of cfChIP vs chromosome-matched downsampled whole-blood (WBC) peaks
2. **Leave-one-out overlap** after dropping a cancer type or a collection site and rebuilding cfChIP consensus peaks

These analyses extend the point-estimate EpiMap overlap in `01_data_preprocessing/` (`epimap_overlap_analysis.sh`).

## Reference epigenomes

Scripts use Roadmap 18-state chromatin-state BED files (`*_18_core_K27ac_dense.bed.gz`) and keep user-specified states, for example:

```text
EnhA1,EnhA2,EnhG1,EnhG2
```

States are grepped from each epigenome, merged with BEDTools, and intersected with cfChIP or WBC peak sets. Epigenome display names come from `label.csv` (EID → NAME).

## Scripts

### `scripts/roadmap_overlap_bootstrapping.sh`

Bootstrap overlap counts against each Roadmap enhancer BED.

**Usage:**
```bash
bash scripts/roadmap_overlap_bootstrapping.sh <patterns> [B] [SEED0] [FRAC]
```

**Parameters:**
- `patterns`: comma-separated ChromHMM states (e.g. `"EnhA1,EnhA2,EnhG1,EnhG2"`)
- `B`: bootstrap replicates (default: 200)
- `SEED0`: base seed (default: 1)
- `FRAC`: within-chromosome sampling fraction for cfChIP (default: 0.8)

**Method:**
1. Build (or reuse) merged enhancer BEDs for the requested states
2. Split cfChIP and WBC consensus peaks by chromosome (drop X, Y, M)
3. For each replicate:
   - **cfChIP**: sample `FRAC` of peaks **with replacement** within each chromosome
   - **WBC**: downsample **without replacement** to the same per-chromosome counts
   - Count unique overlaps (`bedtools intersect -u`) against each reference epigenome
4. Summarize the replicate distribution in R

Reproducible `shuf` streams use `openssl enc -aes-256-ctr` as `--random-source`.

**Output:**
- `work/<states>_combined_results.bootstrap.tsv`
- Long replicate table under `work/tmp/` (`rep`, `EID`, `WB_count`, `CF_count`)

The combined table keeps the `01_data_preprocessing` columns (`EID`, `WB_count`, `CF_count`, `Total_count`, `NAME`), using bootstrap **means** for the count columns, plus `n_boot` and 95% CIs (`WB_sd`, `WB_ci_lo`, `WB_ci_hi`, `CF_sd`, `CF_ci_lo`, `CF_ci_hi`).

### `scripts/summarize_epimap_bootstrap.R`

Summarizes the long bootstrap table. Invoked automatically:

```bash
Rscript scripts/summarize_epimap_bootstrap.R \
  <boot.long.tsv> <total_counts.tsv> <name.tsv> <combined_output.tsv>
```

Place this file in `scripts/` **relative to the working directory** used by `roadmap_overlap_bootstrapping.sh` (the script `cd`s into its analysis folder before calling R).

### `scripts/leave_one_out_overlap.sh`

Leave-one-cancer-type-out and leave-one-collection-site-out overlap.

**Usage:**
```bash
bash scripts/leave_one_out_overlap.sh <patterns> [seed]
```

**Parameters:**
- `patterns`: comma-separated ChromHMM states
- `seed`: WBC downsampling seed (default: 42)

**Method:**
1. Read sample metadata from an Excel sheet (`source`, `antibody`, `study_name`, `cancer_type`, `ctDNA`) and keep H3K27ac samples
2. Match each sample to a peak BED in the cfChIP peak directory
3. Build enhancer BEDs (same as the bootstrap script)
4. For the **full** sample set, and for each leave-out group:
   - Rebuild cfChIP consensus peaks (extend 100 bp, merge, keep intervals overlapping ≥ 5 samples)
   - Downsample WBC consensus peaks per chromosome to match cfChIP counts
   - Count overlaps with each Roadmap enhancer BED
5. Skip a leave-out group if fewer than 10 samples remain
6. Write summary tables and a median cfChIP−WBC difference plot

**Output (under `work/leave_one_out/<states>/`):**
- `leave_one_out_results.tsv`: per-epigenome overlap for each analysis
- `leave_one_out_summary.tsv`: median % overlap and cfChIP−WBC difference
- `leave_one_out_detailed.tsv`: merged EID names and percentages
- `leave_one_out_median_diff.pdf`: bar plot of median differences

## Paths

Both shell scripts currently `cd` into the original cluster working directory and use absolute paths for:

- Roadmap / EpiMap 18-state files
- cfChIP and WBC consensus / per-sample peak BEDs
- genome file for `bedtools slop`
- sample-info Excel workbook
- EpiMap `label.csv`

Update those paths (or copy the scripts next to your data layout) before running. Script logic was not changed.

`roadmap_overlap_bootstrapping.sh` expects `scripts/summarize_epimap_bootstrap.R` inside that working directory.

## Dependencies

- BEDTools (v2.31+ as loaded in the scripts; 2.30+ is sufficient)
- OpenSSL (reproducible `shuf --random-source`)
- R with `data.table` and `ggplot2`
- Python 3 with `openpyxl` (leave-one-out sample-sheet parsing)
- Environment modules as used on the cluster: `module load bedtools/2.31.0`

## Relation to existing modules

| Module | What it estimates |
|--------|-------------------|
| `01_data_preprocessing/epimap_overlap_analysis.sh` | Single overlap point estimate vs Roadmap states |
| `11_roadmap_bootstrap_loo` bootstrap | Uncertainty of those overlap counts |
| `11_roadmap_bootstrap_loo` leave-one-out | Sensitivity to one cancer type or collection site |

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu
