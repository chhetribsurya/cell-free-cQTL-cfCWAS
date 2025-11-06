# Circulating chromatin reveals effects of disease-associated variants on gene regulation

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![BEDTools](https://img.shields.io/badge/BEDTools-2.30%2B-green.svg)](https://bedtools.readthedocs.io/)
[![LDSC](https://img.shields.io/badge/LDSC-1.0.1-orange.svg)](https://github.com/bulik/ldsc)

## Description

This repository contains computational scripts and analysis pipelines for identifying and characterizing genetic determinants of circulating chromatin from large-scale cell-free chromatin profiles. This work presents a comprehensive map of cell-free chromatin QTLs (cfcQTLs) derived from 751 samples using cfChIP-seq targeting H3K27ac (active enhancers and promoters) and H3K4me3 (active promoters).

The repository implements analyses for detecting allelically imbalanced chromatin activity at variant sites, identifying genetic variants that exhibit significant allelic bias in chromatin accessibility and histone modification patterns. Allelic imbalance detection is performed using beta-binomial statistical tests to assess deviations from expected allelic ratios, with multiple testing correction applied via Q-value estimation. This approach identifies allelic imbalance peaks—chromatin regions with significant allelic bias at variant sites—providing insights into cis-regulatory mechanisms and allele-specific chromatin activity.

The repository includes enrichment testing, developmental specificity assessment, heritability estimation using stratified linkage disequilibrium score regression (S-LDSC), cell-free cistrome-wide association studies (cfCWAS) for mapping chromatin-trait associations, and allelic imbalance analysis for identifying allele-specific chromatin activity.

## Overview

This repository is organized into ten main analysis modules:

1. **Setup** (`00_setup/`): Installation instructions and environment configuration
2. **Data Preprocessing** (`01_data_preprocessing/`): Scripts for calculating tissue-specific scores and preparing annotation files
3. **Enrichment Analysis** (`02_enrichment_analysis/`): Scripts for performing cQTL enrichment analyses with permutation testing
4. **LDSC Analysis** (`03_ldsc_analysis/`): Scripts for running stratified LDSC analysis
5. **Tissue Specificity Analysis** (`04_tissue_specificity/`): Scripts for analyzing and visualizing tissue-specific characteristics and comparisons
6. **Genomic Track Visualization** (`05_genomic_tracks/`): Scripts for creating multi-track genomic visualizations
7. **Allelic Imbalance Visualization** (`06_allelic_imbalance/`): Scripts for visualizing allelic imbalance and Q-value analyses
8. **Enrichment Visualization** (`07_enrichment_visualization/`): Scripts for visualizing enrichment analysis results
9. **Manhattan Plot Visualization** (`08_manhattan_plots/`): Scripts for generating Manhattan plots from cfCWAS results
10. **cfCWAS Workflow** (`09_cfcwas_workflow/`): Complete Snakemake workflow for cell-free cistrome-wide association studies

## Repository Structure

```
.
├── 00_setup/                         # Setup and installation instructions
├── 01_data_preprocessing/            # Data preprocessing and annotation scripts
├── 02_enrichment_analysis/           # cQTL enrichment analysis scripts
├── 03_ldsc_analysis/                 # Stratified LDSC analysis scripts
├── 04_tissue_specificity/            # Tissue-specificity analysis and comparisons
├── 05_genomic_tracks/                # Multi-track genomic visualizations
├── 06_allelic_imbalance/             # Allelic imbalance and Q-value visualization
├── 07_enrichment_visualization/      # Enrichment analysis visualization
├── 08_manhattan_plots/               # Manhattan plot generation
└── 09_cfcwas_workflow/               # Cell-free cistrome-wide association studies workflow
```

## Requirements

### Software Dependencies

- **Python 3.x** with the following packages:
  - pandas
  - numpy
  - scipy
  - matplotlib
  - seaborn
  - pybedtools
  - tqdm
  - upsetplot
  - matplotlib_venn

- **R** (version 4.0+) with the following packages:
  - ggplot2
  - dplyr
  - tidyr
  - data.table
  - RColorBrewer
  - Gviz
  - GenomicRanges
  - qvalue
  - Hmisc
  - VennDiagram

- **Command-line tools**:
  - BEDTools (v2.30.0 or higher)
  - LDSC (for stratified LDSC analysis)

### Data Requirements

- cQTL summary statistics files (BED format)
- Roadmap Epigenomics Project chromatin state annotations
- Reference genome files (hg19)
- 1000 Genomes Project reference data (for LDSC)

See `00_setup/README.md` for detailed installation and setup instructions.

## Quick Start

### 0. Setup and Installation

Before running analyses, set up your computational environment:

```bash
cd 00_setup
# Follow instructions in README.md to:
# - Install required software
# - Configure environment variables
# - Download reference data
```

See `00_setup/README.md` for complete setup instructions.

### 1. Data Preprocessing

Calculate tissue-specific scores for chromatin peaks:

```bash
cd 01_data_preprocessing
bash calculate_tissue_specific_score.sh <query_peaks.bed> <assay_name>
```

### 2. Enrichment Analysis

Run cQTL enrichment analysis:

```bash
cd 02_enrichment_analysis
python cqtl_enrichment_analysis.py \
    --cqtl_file <cqtl_file.bed> \
    --random_background <random_background.bed> \
    --regulatory_regions <regulatory_regions.bed> \
    --output_dir <output_directory>
```

### 3. LDSC Analysis

Generate LDSC annotations and run stratified LDSC:

```bash
cd 03_ldsc_analysis
bash make_ldsc_annotation.sh <peaks.bed> <annotation_name> <output_dir>
bash run_stratified_ldsc.sh <trait> <model_files> <baseline_model>
```

### 4. Visualization

Generate figures using the R scripts organized by visualization type:

- **Tissue Specificity** (`04_tissue_specificity/`): Tissue-specificity analysis and comparison plots
- **Genomic Tracks** (`05_genomic_tracks/`): Multi-track genomic visualizations
- **Allelic Imbalance** (`06_allelic_imbalance/`): Q-value distributions and allelic fraction analyses
- **Enrichment Visualization** (`07_enrichment_visualization/`): Regulatory and eQTL enrichment plots
- **Manhattan Plots** (`08_manhattan_plots/`): Genome-wide association plots
- **cfCWAS Workflow** (`09_cfcwas_workflow/`): Complete workflow for cell-free cistrome-wide association studies

Each directory contains scripts organized by analysis type rather than figure numbers.

### 5. cfCWAS Workflow

Run the complete cell-free cistrome-wide association studies workflow:

```bash
cd 09_cfcwas_workflow
conda env create -f env/cfcwas_env.yml
conda activate cfcwas
# Update config.yaml with your data paths
sbatch submit.sh
```

See `09_cfcwas_workflow/README.md` for detailed instructions.

## Analysis Workflow

1. **Data Preprocessing**: Calculate tissue-specific chromatin peak scores by overlapping with Roadmap Epigenomics annotations
2. **Allelic Imbalance Detection**: Identify allelically imbalanced chromatin activity at variant sites using beta-binomial tests, detecting chromatin peaks with significant allelic bias and quantifying allele-specific chromatin activity
3. **Enrichment Analysis**: Test for enrichment of cQTLs in developmentally regulated regulatory regions using Fisher's exact tests and permutation testing
4. **LDSC Analysis**: Estimate heritability enrichment using stratified LDSC
5. **cfCWAS Analysis**: Perform cell-free cistrome-wide association studies to identify genetic associations with chromatin features
6. **Visualization**: Generate publication-quality figures summarizing the results

## Statistical Methods

### Enrichment Testing

Enrichment analyses use Fisher's exact tests to compare overlap frequencies between cQTL sets and size-matched random backgrounds.

### Random Background Generation

Random background sets are generated by sampling from the full pool of available cQTLs that passed equivalent significance thresholds. This approach maintains comparable statistical power and preserves genomic distribution characteristics.

### Allelic Imbalance Detection

Allelic imbalance detection identifies genetic variants exhibiting significant allelic bias in chromatin accessibility and histone modification patterns. The analysis uses beta-binomial statistical tests to assess deviations from expected allelic ratios (50:50) at heterozygous variant sites. Multiple testing correction is applied using Q-value estimation (Storey and Tibshirani, 2003). Allelic imbalance peaks are defined as chromatin regions with significant allelic bias (Q-value < 0.05) at variant sites, providing insights into cis-regulatory mechanisms and allele-specific chromatin activity.

### Regulatory Region Classification

Regulatory regions are classified into four mutually exclusive specificity categories:
- Fetal-specific: Active exclusively in fetal tissues
- Stem cell-specific: Active exclusively in stem cell epigenomes
- Developmental-specific: Active in fetal and/or stem cell tissues but absent in adult tissues
- Adult-specific: Active exclusively in adult tissues

## Citation

If you use these scripts in your research, please cite the associated publication.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

## Documentation

Each analysis module contains its own README with specific usage instructions.

