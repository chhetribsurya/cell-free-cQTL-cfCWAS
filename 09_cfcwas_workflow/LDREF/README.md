# LD Reference Data Directory

This directory contains linkage disequilibrium (LD) reference data from the 1000 Genomes Project, used for LD score regression and association analysis.

## Contents

The directory contains LD reference files for European populations (1000G.EUR):
- **1000G.EUR.{1-22}.bed**: BED files for each chromosome
- **1000G.EUR.{1-22}.bim**: BIM files (variant information) for each chromosome
- **1000G.EUR.{1-22}.fam**: FAM files (sample information) for each chromosome
- **1000G.EUR.snplist**: List of SNPs in the reference panel
- **hm3.pos**: HapMap3 SNP positions

## Purpose

LD reference data is essential for:
- Computing LD scores for stratified LDSC analysis
- Performing association tests with proper LD adjustment
- Estimating heritability and enrichment statistics
- Correcting for population structure and LD patterns

## Data Source

These files are derived from the 1000 Genomes Project and FUSION LD reference data:
- Download URL: https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
- Original source: 1000 Genomes Project Phase 3

## Installation

To set up the LD reference data:

```bash
cd LDREF
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar -xjf LDREF.tar.bz2
```

## Usage

The LD reference files are automatically used by:
- FUSION software for weight computation
- LDSC for heritability estimation
- Association testing scripts

## Note

Ensure that the LD reference data matches the genome build (hg19) used in your analysis. The reference data should be placed in this directory before running the workflow.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

