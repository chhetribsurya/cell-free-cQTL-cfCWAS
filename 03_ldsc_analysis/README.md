# Stratified LDSC Analysis

This directory contains scripts for performing stratified linkage disequilibrium score regression (S-LDSC) analysis to estimate heritability enrichment of cQTLs and regulatory regions.

## Scripts

### `make_ldsc_annotation.sh`

Generates LDSC annotation files from BED files for use in stratified LDSC analysis.

**Usage:**
```bash
bash make_ldsc_annotation.sh <peaks.bed> <annotation_name> <output_directory>
```

**Parameters:**
- `peaks.bed`: Input BED file containing genomic regions (peaks)
- `annotation_name`: Name identifier for the annotation
- `output_directory`: Base directory for output files

**Requirements:**
- LDSC software installed
- 1000 Genomes Project reference data (1000G_EUR_Phase3_plink)
- HapMap3 SNP list

**Method:**
1. Sorts and formats the input BED file
2. Adds "chr" prefix to chromosome names if missing
3. Generates annotation files for each chromosome (1-22) using `ldsc/make_annot.py`
4. Computes LD scores for each chromosome using `ldsc/ldsc.py`

**Output:**
- `<annotation_name>.<chr>.annot.gz`: Annotation files for each chromosome
- `<annotation_name>.<chr>.l2.ldscore.gz`: LD score files for each chromosome
- `<annotation_name>.<chr>.l2.M`: Number of SNPs in each annotation

**Notes:**
- The script assumes a specific directory structure with LDSC installation and reference data
- Requires SLURM job scheduler (can be modified for local execution)
- Processes chromosomes 1-22

### `run_stratified_ldsc.sh`

Runs stratified LDSC analysis to estimate heritability enrichment for a given trait.

**Usage:**
```bash
bash run_stratified_ldsc.sh <trait> <model_files> <baseline_model>
```

**Parameters:**
- `trait`: Name of the trait (corresponds to summary statistics file)
- `model_files`: Path to model annotation files (LD score files)
- `baseline_model`: Name of the baseline model to use

**Input Files:**
- Summary statistics file: `work/sumstats/<trait>.sumstats.gz`
- LD scores: Specified in `model_files` parameter
- Weights file: `data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.`
- Frequency file: `data/1000G_Phase3_frq/1000G.EUR.QC.`

**Output:**
- `<trait>_<baseline_model>.results`: LDSC results file
- `<trait>_<baseline_model>.log`: Log file with detailed output

**Method:**
The script runs LDSC with the following options:
- `--h2`: Estimates heritability using summary statistics
- `--w-ld-chr`: LD weights for regression
- `--ref-ld-chr`: Reference LD scores (model annotations)
- `--overlap-annot`: Accounts for annotation overlap
- `--frqfile-chr`: Allele frequency files
- `--print-coefficients`: Prints coefficients for each annotation

**Output Interpretation:**
The results file contains:
- Heritability estimates
- Enrichment statistics for each annotation
- Standard errors and p-values
- Coefficient estimates

## LDSC Background

Stratified LDSC is a method for partitioning heritability across genomic annotations. It estimates the contribution of different genomic regions to trait heritability by:
1. Computing LD scores for each annotation
2. Regressing summary statistics on LD scores
3. Estimating enrichment as the ratio of per-SNP heritability in the annotation to the genome-wide average

## Requirements

- LDSC software (available from https://github.com/bulik/ldsc)
- Python 2.7 (for LDSC)
- 1000 Genomes Project Phase 3 reference data
- HapMap3 SNP list
- BEDTools for preprocessing

## Setup

1. Install LDSC following the official documentation
2. Download and prepare 1000 Genomes Project reference data
3. Prepare summary statistics files in LDSC format
4. Generate annotation files using `make_ldsc_annotation.sh`
5. Run stratified LDSC using `run_stratified_ldsc.sh`

## Notes

- LDSC requires Python 2.7 (legacy software)
- Summary statistics must be in LDSC format
- All files must use the same genome build (typically hg19/GRCh37)
- The scripts are configured for SLURM job scheduler but can be adapted for local execution
- Processing is done per chromosome for efficiency

## Citation

If using LDSC, please cite:
- Finucane HK, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nature Genetics. 2015.

