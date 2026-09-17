# Covariate-Matched Random Backgrounds

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![nullranges](https://img.shields.io/badge/nullranges-matchRanges-green.svg)](https://bioconductor.org/packages/nullranges)
[![BEDTools](https://img.shields.io/badge/BEDTools-2.30%2B-green.svg)](https://bedtools.readthedocs.io/)
[![Bash](https://img.shields.io/badge/Bash-4.0%2B-lightgrey.svg)](https://www.gnu.org/software/bash/)

This directory generates genome-wide **covariate-matched null ranges** for cQTL / cfChIP enrichment testing using Bioconductor `nullranges::matchRanges()`, then runs SLURM-array enrichment against target annotations (eQTLs, asaQTLs, or other BED sets).

It extends the permutation enrichment in `02_enrichment_analysis/` and `09_cfcwas_workflow/scripts/` by replacing length-only `bedtools shuffle` backgrounds with a matched null set.

## What was used: `matchRanges`

Matching follows the `nullranges` `matchRanges()` workflow described in the Bioconductor vignette:

- Vignette: [Introduction to matchRanges](https://www.bioconductor.org/packages/release/bioc/vignettes/nullranges/inst/doc/matchRanges.html)
- Package: [nullranges](https://bioconductor.org/packages/nullranges)
- Method paper: Davis et al. (2023), *Bioinformatics* ([doi:10.1093/bioinformatics/btad197](https://doi.org/10.1093/bioinformatics/btad197))

`matchRanges()` selects a `matched` subset of a `pool` that has a similar covariate distribution to a `focal` set, using propensity scores. Supported methods are nearest-neighbor, rejection sampling, and stratified sampling (with or without replacement).

This module uses **stratified sampling without replacement**, the recommended alternative because nearest-neighbor matching without replacement is not implemented in `matchRanges()`.

### Covariates matched here

| Covariate | How it is computed | Role |
|-----------|--------------------|------|
| Chromosome | `bedtools shuffle -chrom` plus a chromosome factor in `matchRanges` | Preserve chromosome distribution |
| Length | Shuffle of length-matched foreground intervals | Same interval size as focal peaks |
| GC content | `bedtools nuc` | Sequence composition |
| Distance to TSS | `bedtools closest` to RefSeq TSS, then `log1p` | Gene proximity (heavy-tailed) |
| Mappability | `bigWigAverageOverBed` mean0 | Mapping uniqueness |
| Accessibility | `bigWigAverageOverBed` mean0 | Open-chromatin / alignability track |
| Interval length | `end - start` | Peak width |
| N-base fraction | `bedtools nuc`; candidates with `pctN >= 0.05` are dropped | Avoid assembly gaps |

## Scripts

### `run_code.sh`

Driver that writes a target BED list and submits `scripts/enrichment_array.sh` as a SLURM array.

**Usage:**
```bash
# Edit paths in run_code.sh, then:
bash run_code.sh
```

**Configuration (edit at the top of the script):**
- `FG`: foreground BED (e.g. significant cfChIP peaks / cQTLs)
- `BG_UNIV`: background universe BED (consensus peaks) for the existing permutation null
- `TARGET_DIR`: directory of target annotation BEDs
- `GENOME_SIZES`, `FASTA`, `TSS`, `BLACKLIST`, `MAPBW`, `ACCBW`: hg19 resources
- `NBATCHES`: permutation / matched-random batches (default: 50)
- `CAND_MULT`: candidate-pool multiplier for `matchRanges` (default: 50)

The example in the script uses conda environment `enrich-matched-bg` and hg19 resources. Update cluster paths before running.

### `scripts/make_matched_random_bg.sh`

Builds a chromosome- and length-matched candidate pool with `bedtools shuffle`, computes covariates, and calls `make_matched_random_bg.R`.

**Usage:**
```bash
bash scripts/make_matched_random_bg.sh \
  <fg.bed> \
  <genome.sizes.nochr> \
  <fasta.fa> \
  <tss.bed.nochr> \
  <blacklist.bed.nochr> \
  <mappability.bw> \
  <accessibility.bw> \
  <seed> \
  <multiplier> \
  <out.bed>
```

**Output:**
- `out.bed`: 3-column BED without a `chr` prefix (`chr start end`)

**Method:**
1. Strip `chr` prefixes and add unique foreground IDs
2. Replicate foreground intervals `multiplier` times and `bedtools shuffle -chrom -excl blacklist`
3. Compute GC, pctN, TSS distance, mappability, accessibility, and length
4. Drop intervals with N-base fraction ≥ 0.05
5. Run `matchRanges()` and write the matched BED

### `scripts/make_matched_random_bg.R`

R implementation of propensity-score matching.

**Usage:**
```bash
Rscript scripts/make_matched_random_bg.R <fg.features.tsv> <cand.features.tsv> <out.bed>
```

**Input feature tables** (tab-separated, no header):

```
chr  start  end  name  gc  pctN  distTSS  map  acc  len
```

**`matchRanges` call used:**
```r
matchRanges(
  focal   = fg_gr,
  pool    = cand_gr,
  covar   = ~ gc + logDistTSS + map + acc + len + chrF,
  method  = "stratified",
  replace = FALSE
)
```

### `scripts/enrichment_array.sh`

SLURM array job. For each target BED it:

1. Computes foreground enrichment (`scripts/enrich.sh`)
2. Runs universe-background permutations (`scripts/enrich.permute.sh`)
3. Generates one matched-random background per batch and repeats permutations
4. Writes empirical p-values, null mean/SD/SE, and relative enrichment with confidence intervals for both nulls

**Usage (via `run_code.sh` or sbatch):**
```bash
sbatch --array=1-${NTARGETS} \
  --export=ALL,FG=...,BG_UNIV=...,OUTDIR=...,TARGET_LIST=..., \
           GENOME_SIZES=...,FASTA=...,TSS=...,BLACKLIST=..., \
           MAPBW=...,ACCBW=...,CAND_MULT=50,TMPROOT=/tmp,NBATCHES=50 \
  scripts/enrichment_array.sh
```

**Output (per target):**
- `work/output/<OUTDIR>/<group>/<name>/fg/fg.enrich.txt`
- `perm/group.*.bg.enrich`: universe-background null
- `rand/bg.matched.seed*.bed` and `rand/group.*.bg.enrich`: matched-random null
- `enrichment.stats.txt` and `enrichment.stats.with_fg.txt`

### Helper scripts (same enrichment engine as `09_cfcwas_workflow/scripts/`)

| File | Role |
|------|------|
| `scripts/enrich.sh` | Foreground overlap enrichment (bp overlap × 1000 / total bp) |
| `scripts/enrich.permute.sh` | Random draws from a background BED, same enrichment metric |
| `scripts/sum.awk` | Sum a numeric column |
| `scripts/mean_sd_se.awk` | Mean, SD, SE, and N of a null distribution |
| `scripts/quantile.awk` | Empirical quantile of a sorted stream (R type-7) |

`enrich.sh` and `enrich.permute.sh` are copies of the existing cfCWAS enrichment helpers so this module is self-contained. They were not modified.

## Expected directory layout when running

`run_code.sh` and `enrichment_array.sh` resolve helpers as `scripts/...` from the working directory:

```
10_matched_random_background/
├── run_code.sh
├── scripts/
│   ├── enrichment_array.sh
│   ├── make_matched_random_bg.sh
│   ├── make_matched_random_bg.R
│   ├── enrich.sh
│   ├── enrich.permute.sh
│   ├── sum.awk
│   ├── mean_sd_se.awk
│   └── quantile.awk
├── data/                         # not in repo; supply locally
│   ├── foreground/
│   ├── background/
│   ├── target/
│   └── resources/
└── work/
```

## Dependencies

- BEDTools (v2.30+)
- `bigWigAverageOverBed` (UCSC kent utilities)
- R ≥ 4.0 with Bioconductor packages:
  - `nullranges`
  - `GenomicRanges`
  - `S4Vectors`
- SLURM (for the array driver)
- conda environment `enrich-matched-bg` (as named in `run_code.sh`)

Install `nullranges`:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("nullranges")
```

## Notes

- Cluster paths in `run_code.sh` (reference FASTA, chrom.sizes, SLURM account email) are from the original analysis environment. Update them for your system; do not change script logic.
- Chromosome names in genome, TSS, and blacklist files should be **without** a `chr` prefix (`*.nochr`), matching `make_matched_random_bg.sh`.
- `enrichment_array.sh` passes an extra `TMPROOT` argument to `make_matched_random_bg.sh`; the matcher uses `${TMPDIR:-/tmp}` internally and ignores that extra argument.
- Random seed is fixed at 1 inside `make_matched_random_bg.R` after candidate generation; shuffle seeds still vary by batch via the shell `seed` argument.

## References

1. Davis ES, Mu W, Lee S, Dozmorov MG, Love MI, Phanstiel DH. matchRanges: generating null hypothesis genomic ranges via covariate-matched sampling. *Bioinformatics* (2023). https://doi.org/10.1093/bioinformatics/btad197
2. `nullranges` `matchRanges` vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/nullranges/inst/doc/matchRanges.html

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu
