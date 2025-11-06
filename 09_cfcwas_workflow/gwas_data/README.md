# GWAS Data Directory

This directory contains GWAS (Genome-Wide Association Study) summary statistics files used for association testing in the cfCWAS workflow.

## Contents

GWAS summary statistics files should be placed in this directory with the naming convention `*.sumstats.gz` (gzipped tab-delimited files).

## File Format

Each GWAS summary statistics file should have the following format:

```
SNP    A1    A2    N      Z
rs123  A     G     100000 2.5
rs456  T     C     100000 -1.8
```

Where:
- **SNP**: SNP identifier (rsID)
- **A1**: Effect allele
- **A2**: Alternative allele
- **N**: Sample size
- **Z**: Z-score (effect size / standard error)

## Purpose

GWAS summary statistics are used to:
- Test associations between predicted chromatin activity and traits
- Identify chromatin features genetically associated with disease/trait risk
- Perform cistrome-wide association studies (cfCWAS)

## Usage

1. Place your GWAS summary statistics files in this directory
2. Name files with the pattern `{trait_name}.sumstats.gz`
3. Update the `gwas_traits` list in `config.yaml` to include your trait names
4. The workflow will automatically process these files during association testing

## Example Files

The directory may contain example GWAS files from various studies:
- Disease traits (e.g., prostate cancer, type 2 diabetes)
- Quantitative traits (e.g., testosterone levels)
- UK Biobank phenotypes

## Note

- Files must be gzipped (`.gz` extension)
- Ensure SNP identifiers match the reference panel used (typically rsIDs)
- The workflow expects standard GWAS summary statistics format
- Large files may take time to process

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

