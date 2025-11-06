# GitHub Repository

Computational pipelines for cell-free chromatin QTL (cfcQTL) mapping, allelic imbalance detection, and cistrome-wide association studies (cfCWAS) from cfChIP-seq data.

## Description

This repository provides computational scripts and analysis pipelines for identifying genetic determinants of circulating chromatin from large-scale cell-free chromatin profiles. The repository implements:

- **Cell-free chromatin QTL (cfcQTL) mapping** from cfChIP-seq data (H3K27ac and H3K4me3)
- **Allelic imbalance detection** using beta-binomial statistical tests to identify variant-specific chromatin activity
- **Enrichment analysis** for developmental specificity and regulatory region classification
- **Stratified LDSC** for heritability estimation
- **Cell-free cistrome-wide association studies (cfCWAS)** for mapping chromatin-trait associations

The analysis pipeline includes data preprocessing, enrichment testing, visualization, and a complete Snakemake workflow for cfCWAS. This liquid biopsy approach enables mapping of genetic effects on regulatory element activity in tissues typically inaccessible through conventional sampling.

**Tools and Softwares**: Python, R, Snakemake, BEDTools, LDSC

## Alternative Shorter Version

Computational pipelines for cell-free chromatin QTL (cfcQTL) mapping and cistrome-wide association studies (cfCWAS) from cfChIP-seq data. Includes allelic imbalance detection, enrichment analysis, heritability estimation, and visualization tools. Implemented in Python, R, and Snakemake.

