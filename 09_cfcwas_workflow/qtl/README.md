# QTL Data Directory

This directory contains quantitative trait locus (QTL) data files in BED format, representing chromatin QTLs mapped from various tissues and cell types.

## Contents

The directory contains QTL files for multiple tissues and cell types, each named as `{Tissue_Name}.qtl.hg19.bed`. These files contain:
- Chromosome coordinates (chr, start, end)
- QTL effect sizes and statistics
- SNP identifiers and genomic positions

## File Format

Each QTL file is in BED format with additional columns containing QTL statistics:
- Standard BED columns: chromosome, start, end
- Additional columns: QTL effect sizes, p-values, SNP information

## Purpose

These QTL files are used to:
1. Compute prediction weights for chromatin features
2. Model chromatin activity as a function of nearby SNP genotypes
3. Predict chromatin activity for individuals in GWAS studies

## Usage

The QTL files are automatically loaded by the Snakemake workflow during the heritability analysis and weight computation phases. Users should ensure that QTL files match the tissue types and chromatin marks being analyzed in their cfCWAS study.

## Note

These QTL files are reference data from the original CWAS workflow. For cell-free chromatin analysis, users may need to generate their own QTL files from cfChIP-seq data or use appropriate reference QTLs that match their experimental design.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

