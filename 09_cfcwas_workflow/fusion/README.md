# FUSION Integration Directory

This directory contains scripts and reference files for integrating FUSION (Functional Summary-based Imputation) software into the cfCWAS workflow.

## Contents

- **FUSION.assoc_test.R**: R script for performing association testing using FUSION
- **FUSION.post_process.R**: R script for post-processing FUSION results
- **glist-hg19**: Gene list file for hg19 genome build

## Purpose

FUSION is used for transcriptome-wide association studies (TWAS) and provides methods for:
- Computing prediction weights for chromatin features
- Performing association tests between predicted chromatin activity and traits
- Post-processing and summarizing association results

## Dependencies

- FUSION software: https://github.com/gusevlab/fusion_twas
- R packages: data.table, dplyr, and other dependencies specified in the workflow

## Usage

These scripts are called by the Snakemake workflow (`cfcwas.snakefile`) during the association testing phase of the analysis.

## Reference

FUSION software and documentation: https://github.com/gusevlab/fusion_twas

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

