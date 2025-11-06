# Manhattan Plot Visualization

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![ggplot2](https://img.shields.io/badge/ggplot2-3.3%2B-green.svg)](https://ggplot2.tidyverse.org/)
[![data.table](https://img.shields.io/badge/data.table-1.14%2B-blue.svg)](https://rdatatable.gitlab.io/data.table/)
[![dplyr](https://img.shields.io/badge/dplyr-1.0%2B-blue.svg)](https://dplyr.tidyverse.org/)

This directory contains R scripts for generating Manhattan plots from cQTL and CWAS (Cell-free chromatin Wide Association Study) analysis results. These plots visualize genome-wide association signals across chromosomes for multiple traits.

## Scripts

### `manhattan_plot_cancer_traits.R`

Creates Manhattan plots for cancer-related traits from CWAS analysis results.

**Input Files:**
- CWAS results files from `Results/CWAS/Results/cfChIP_H3K4me3/Cancer/` directory
- Dummy SNPs file for preserving full chromosome spans (`Results/dummysnps.rds`)

**Output:**
- Manhattan plot showing association signals across chromosomes
- Multiple traits displayed in a single or multi-panel plot
- PDF figure saved to manuscript figures directory

**Key Features:**
- Displays TWAS (Transcriptome-Wide Association Study) p-values
- Orders traits by number of significant signals
- Color-coded by chromosome
- Includes significance thresholds
- Shows cumulative genomic positions for proper chromosome spacing

### `manhattan_plot_non_cancer_traits.R`

Creates Manhattan plots for non-cancer traits from CWAS analysis results.

**Input Files:**
- CWAS results files from `Results/CWAS/Results/cfChIP_H3K4me3/Non-Cancer/` directory
- Dummy SNPs file

**Output:**
- Manhattan plot for non-cancer traits
- PDF figure saved to manuscript figures directory

**Key Features:**
- Similar to cancer traits but focused on non-cancer phenotypes
- Filters out specific traits (e.g., Psoriasis) as needed
- Same visualization approach as cancer traits

### `manhattan_plot_cwas_results.R`

Generates Manhattan plots from general CWAS analysis results.

**Input Files:**
- CWAS results files
- Dummy SNPs file

**Output:**
- Manhattan plot showing CWAS association signals
- PDF figure saved to manuscript figures directory

**Key Features:**
- Displays genome-wide association results
- Chromosome-wise visualization
- Significance thresholds

### `manhattan_plot_cwas_results_k27.R`

Creates Manhattan plots specifically for H3K27ac CWAS results.

**Input Files:**
- H3K27ac-specific CWAS results files
- Dummy SNPs file

**Output:**
- Manhattan plot for H3K27ac associations
- PDF figure saved to manuscript figures directory

**Key Features:**
- Histone modification-specific visualization
- H3K27ac association patterns

## Dependencies

- R packages: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `RColorBrewer`

## Installation

Install required packages:

```r
install.packages(c("data.table", "dplyr", "tidyr", "ggplot2", "RColorBrewer"))
```

## Usage

Each script is self-contained and can be run independently:

```bash
Rscript manhattan_plot_cancer_traits.R
Rscript manhattan_plot_non_cancer_traits.R
Rscript manhattan_plot_cwas_results.R
Rscript manhattan_plot_cwas_results_k27.R
```

**Note:** Update file paths in each script to match your local directory structure. Look for `setwd()` calls and file paths at the beginning of each script.

## Manhattan Plot Features

### Chromosome Visualization

- X-axis: Cumulative genomic position across all chromosomes
- Y-axis: -log10(p-value) for association signals
- Colors: Alternating colors for different chromosomes
- Labels: Chromosome numbers at chromosome centers

### Significance Thresholds

- Horizontal lines indicate significance thresholds
- Trait-specific thresholds based on significant signal counts
- Bonferroni correction or FDR thresholds as appropriate

### Data Processing

1. **Dummy SNPs**: Used to preserve full chromosome spans and ensure proper spacing
2. **Trait Ordering**: Traits are ordered by number of significant signals (ascending)
3. **Position Calculation**: Midpoint positions calculated from P0 and P1 coordinates
4. **Signal Identification**: Significant signals identified based on trait-specific thresholds

## Statistical Methods

### P-value Calculation

- TWAS p-values are used for association signals
- Transformed to -log10 scale for visualization
- Significance determined by trait-specific thresholds

### Multiple Testing Correction

- Bonferroni correction may be applied
- FDR thresholds used where appropriate
- Trait-specific significance levels

## Output Format

All scripts generate PDF figures with:
- Consistent color schemes across chromosomes
- High-resolution vector graphics
- Standardized figure dimensions
- Proper chromosome labeling and spacing

## Customization

### Color Schemes

Modify chromosome colors in the script:
```r
chromosome_colors <- c("#1f78b4", "#33a02c", ...)  # Define color palette
```

### Significance Thresholds

Adjust significance lines:
```r
sigline <- max(model.sig$TWAS.P, na.rm = TRUE)  # Trait-specific threshold
```

### Figure Dimensions

Modify output dimensions:
```r
ggsave(..., width = 200, height = 150, units = "mm")
```

## Notes

- All coordinates are in hg19/GRCh37
- Dummy SNPs are required for proper chromosome visualization
- Scripts handle missing data gracefully
- Multiple traits can be displayed in single or multi-panel layouts

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

