# Allelic Imbalance Analysis Visualization

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![qvalue](https://img.shields.io/badge/qvalue-2.26%2B-green.svg)](https://bioconductor.org/packages/qvalue/)
[![ggplot2](https://img.shields.io/badge/ggplot2-3.3%2B-green.svg)](https://ggplot2.tidyverse.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14%2B-green.svg)](https://bioconductor.org/)

This directory contains R scripts for visualizing allelic imbalance (AI) analysis results, including Q-value distributions, cQTL peak proportions, allelic fraction correlations, and MPRA overlap visualizations.

## Scripts

### `qvalue_distribution_histogram.R`

Generates histograms showing the distribution of log10-transformed read-count ratios between cfChIP and WBC for significantly imbalanced peaks.

**Input Files:**
- Combined results file from allelic imbalance analysis (`Results/cf_wb_AI/cf_peaks/results.all.txt`)

**Output:**
- Histogram of log10(cfChIP/WBC) ratios for cfChIP-significant peaks
- Histogram for all tested peaks
- PDF figures saved to `Manuscripts/Figures/ext_2a.pdf` and `Figures/cf_wb_AI/p.count.dist.all.pdf`

**Key Features:**
- Filters peaks by Q-value thresholds (BBINOM.C1.Q < 0.05, BBINOM.C0.Q >= 0.05)
- Shows mean ratio as vertical line
- Displays zero line reference

### `qvalue_bin_statistical_analysis.R`

Performs statistical analysis of Q-value distributions across bins and generates boxplots with statistical testing.

**Input Files:**
- Combined results file from allelic imbalance analysis

**Output:**
- Boxplot showing Q-value distributions across bins
- Statistical test results (ANOVA and non-parametric tests)
- PDF figure with statistical annotations

**Key Features:**
- Categorizes results into bins based on statistical significance
- Performs ANOVA and non-parametric statistical tests
- Includes subsampling for large non-significant bins
- Displays p-values and statistical annotations

### `cqtl_peak_proportion_pie.R`

Creates pie charts showing the proportion of cQTL peaks relative to total peaks for H3K27ac and H3K4me3 assays.

**Input Files:**
- Peak count data (hardcoded in script)

**Output:**
- Pie charts for H3K27ac and H3K4me3 showing:
  - Total peaks
  - Imbalanced/cQTL peaks
- PDF figures saved to `Manuscripts/Figures/2_b_pie_H3K27ac.pdf` and `2_b_pie_H3K4me3.pdf`

**Key Features:**
- Shows proportions of cQTL peaks
- Color-coded by peak type
- Displays percentages

### `allelic_fraction_correlation_scatter.R`

Creates scatter plots showing the correlation of allelic fractions between cfChIP and WBC ChIP-seq data.

**Input Files:**
- cfChIP allelic imbalance summary statistics (`Results/AI_cQTLs/cfChIP.k27.AI.all.sumstats.txt`)
- WBC allelic imbalance summary statistics (`Results/AI_cQTLs/wb.k27.AI.all.sumstats.txt`)

**Output:**
- Scatter plot of allelic fractions (cfChIP vs WBC)
- Correlation coefficient and p-value
- PDF figure saved to `Manuscripts/Figures/2_e.pdf`

**Key Features:**
- Shows Pearson correlation
- Color-codes points by allelic fraction category:
  - Both < 0.5 (red)
  - Both > 0.5 (blue)
  - Opposite direction (grey)
- Includes regression line and correlation statistics

### `mpra_overlap_visualization.R`

Visualizes overlap between cQTLs and MPRA (Massively Parallel Reporter Assay) data, showing allelic activity for specific SNPs.

**Input Files:**
- MPRA data file (`Data/MPRA/df_allelic_activity_cancer_SNV.xlsx`)
- cQTL results file (`Results/AI_cQTLs/cfChIP.k4.combined.sig.txt`)
- Sample information files

**Output:**
- Comparison plots showing allelic activity for overlapping SNPs
- PDF figure for specific SNP visualization

**Key Features:**
- Identifies cQTLs that overlap with MPRA-validated SNPs
- Shows allelic activity patterns
- Visualizes specific loci of interest

## Dependencies

- R packages: `qvalue`, `ggplot2`, `tidyr`, `dplyr`, `Hmisc`, `VennDiagram`, `readr`, `data.table`, `ggpubr`, `rstatix`, `forcats`, `readxl`, `scales`, `patchwork`

## Installation

Install required packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "data.table", 
                   "qvalue", "Hmisc", "VennDiagram", "readr", 
                   "ggpubr", "rstatix", "forcats", "readxl", 
                   "scales", "patchwork"))

# Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
```

## Usage

Each script is self-contained and can be run independently:

```bash
Rscript qvalue_distribution_histogram.R
Rscript qvalue_bin_statistical_analysis.R
Rscript cqtl_peak_proportion_pie.R
Rscript allelic_fraction_correlation_scatter.R
Rscript mpra_overlap_visualization.R
```

**Note:** Update file paths in each script to match your local directory structure. Look for `setwd()` calls and file paths at the beginning of each script.

## Statistical Methods

### Q-value Calculation

Q-values are calculated using the `qvalue` package, which implements the method of Storey and Tibshirani (2003) for estimating false discovery rates.

### Correlation Analysis

Pearson correlation is used to assess the relationship between allelic fractions in cfChIP and WBC datasets. P-values are calculated using standard correlation tests.

### Bin Analysis

Statistical tests (ANOVA and non-parametric alternatives) are performed to compare Q-value distributions across different bins of significance.

## Output Format

All scripts generate publication-ready PDF figures with:
- Helvetica font family
- Consistent color schemes
- High-resolution vector graphics
- Statistical annotations where appropriate

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

