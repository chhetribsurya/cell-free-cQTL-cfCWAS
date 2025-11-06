# Heritability Analysis Directory

This directory contains scripts for computing prediction weights and performing heritability analysis using FUSION.

## Contents

- **FUSION.compute_weights.R**: R script for computing prediction weights for chromatin features
- **compute_weights.sh**: Shell script wrapper for weight computation
- **reformat_covar.R**: R script for reformatting covariate files

## Purpose

The heritability analysis (hsq) step computes prediction weights that model chromatin activity as a function of nearby SNP genotypes. These weights are then used to predict chromatin activity for individuals in GWAS studies.

## Workflow

1. **Weight Computation**: Uses FUSION to compute prediction weights from QTL data
2. **Covariate Reformating**: Prepares covariate files for downstream analysis
3. **Integration**: Weights are used in the association testing phase

## Dependencies

- FUSION software
- R packages: data.table, dplyr
- QTL data files (from `qtl/` directory)

## Usage

These scripts are automatically called by the Snakemake workflow during the heritability analysis phase.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

