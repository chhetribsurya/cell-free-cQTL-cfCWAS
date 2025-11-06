# Enrichment Analysis Visualization

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![ggplot2](https://img.shields.io/badge/ggplot2-3.3%2B-green.svg)](https://ggplot2.tidyverse.org/)
[![ComplexHeatmap](https://img.shields.io/badge/ComplexHeatmap-2.10%2B-blue.svg)](https://bioconductor.org/packages/ComplexHeatmap/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14%2B-green.svg)](https://bioconductor.org/)

This directory contains R scripts for visualizing enrichment analysis results, including regulatory region enrichments, GTEx eQTL enrichments, and comparative enrichment analyses across different datasets and histone modifications.

## Scripts

### `regulatory_enrichment_barplots.R`

Creates bar plots showing fold enrichment of cQTLs in different regulatory region categories (fetal-specific, stem cell-specific, developmental-specific, adult-specific).

**Input Files:**
- Detailed statistics file from enrichment analysis (`detailed_statistics_full_random.csv`)

**Output:**
- Bar plots showing fold enrichment by category
- Odds ratio plots with confidence intervals
- P-value plots
- PDF figures saved to enrichment plots directory

**Key Features:**
- Color-coded by regulatory category
- Shows fold enrichment, odds ratios, and statistical significance
- Handles multiple regulatory element types (active, enhancer, promoter, bivalent)

### `regulatory_enrichment_barplots_k4.R` and `regulatory_enrichment_barplots_k27.R`

Variant scripts for H3K4me3 and H3K27ac-specific enrichment visualizations.

**Key Features:**
- Assay-specific enrichment patterns
- Histone modification-specific regulatory associations

### `regulatory_enrichment_comparison.R`

Compares enrichment patterns across different conditions or datasets.

**Key Features:**
- Side-by-side comparisons
- Statistical comparisons between conditions

### `regulatory_enrichment_comparison_wbc.R` and `regulatory_enrichment_comparison_combined.R`

Comparison scripts for WBC-specific and combined analyses.

**Key Features:**
- cfChIP vs WBC comparisons
- Combined visualization of multiple datasets

### `regulatory_enrichment_summary.R`

Generates summary visualizations of regulatory enrichment results.

**Key Features:**
- Comprehensive overview of enrichment patterns
- Summary statistics and key findings

### `epimap_overlap_comparison.R`

Creates horizontal bar plots comparing EpiMap tissue overlap percentages between cfChIP and WBC ChIP data.

**Input Files:**
- EpiMap overlap grouped output ratios file

**Output:**
- Horizontal bar plot showing difference in overlap percentages
- PDF figure saved to `Manuscripts/Figures/3_a.pdf`

**Key Features:**
- Compares cfChIP and WBC ChIP overlap with Roadmap Epigenomics tissues
- Color-coded by direction of difference
- Shows tissue-specific patterns

### `regulatory_enrichment_heatmap.R`

Creates heatmaps showing enrichment patterns across different regulatory categories and conditions.

**Key Features:**
- Multi-dimensional visualization of enrichment
- Color-coded by enrichment strength

### `gtex_eqtl_enrichment_heatmap.R`

Generates heatmaps showing enrichment of cQTLs in GTEx eQTL data across different tissues and conditions.

**Input Files:**
- GTEx eQTL enrichment summary files

**Output:**
- Heatmap showing enrichment across tissues
- PDF figure saved to `Manuscripts/Figures/3_e.pdf`

**Key Features:**
- Tissue-specific enrichment patterns
- Comparison across different cQTL types (H3K4me3, H3K27ac, WBC)
- Color-coded by enrichment strength
- Grouped by tissue categories (Adipose, Artery, Brain, etc.)

### `gtex_eqtl_enrichment_heatmap_k27.R` and `gtex_eqtl_enrichment_heatmap_k4.R`

Variant scripts for H3K27ac and H3K4me3-specific GTEx eQTL enrichment heatmaps.

**Key Features:**
- Histone modification-specific patterns
- Tissue-specific associations

### `enrichment_summary_plots.R`

Generates comprehensive summary plots for cQTL and cfPeak enrichment analyses across multiple regulatory element types.

**Input Files:**
- Enrichment analysis results for multiple regulatory element categories

**Output:**
- Summary plots for active, enhancer, promoter, and bivalent regulatory elements
- Multiple plot types (fold enrichment, odds ratios, p-values)

**Key Features:**
- Comprehensive visualization across regulatory categories
- Multiple statistical metrics
- Publication-ready formatting

## Dependencies

- R packages: `ggplot2`, `dplyr`, `tidyr`, `data.table`, `RColorBrewer`, `ComplexHeatmap`, `circlize`, `grid`, `readr`, `stringr`

## Installation

Install required packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "data.table", 
                   "RColorBrewer", "readr", "stringr", "grid"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))
```

## Usage

Each script is self-contained and can be run independently:

```bash
Rscript regulatory_enrichment_barplots.R
Rscript epimap_overlap_comparison.R
Rscript gtex_eqtl_enrichment_heatmap.R
Rscript enrichment_summary_plots.R
```

**Note:** Update file paths in each script to match your local directory structure. Look for `setwd()` calls and file paths at the beginning of each script.

## Regulatory Categories

Scripts visualize enrichment across four mutually exclusive regulatory categories:

1. **Fetal-specific**: Active exclusively in fetal tissues
2. **Stem cell-specific**: Active exclusively in stem cell epigenomes
3. **Developmental-specific**: Active in fetal and/or stem cell tissues but absent in adult tissues
4. **Adult-specific**: Active exclusively in adult tissues

## Color Schemes

Standard color mappings used across scripts:
- Fetal-specific: `#F1C40F` (Yellow)
- Developmental-specific: `#E74C3C` (Red)
- Stem cell-specific: `#3498DB` (Blue)
- Adult-specific: `#2ECC71` (Green)

## Output Format

All scripts generate publication-ready PDF figures with:
- Helvetica font family
- Consistent color schemes
- High-resolution vector graphics
- Standardized figure dimensions

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

