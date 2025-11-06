# Sample Metadata Visualization

This directory contains R scripts for visualizing sample composition, metadata, and tissue-specific characteristics of the cQTL dataset.

## Scripts

### `sample_composition_by_cancertype.R`

Generates bar plots showing the distribution of samples across cancer types for different histone modification assays (H3K4me3 and H3K27ac).

**Input Files:**
- Sample metadata Excel file
- Sample information files for H3K4me3 and H3K27ac assays

**Output:**
- Bar plot showing sample counts by cancer type and assay type
- PDF figure saved to `Manuscripts/Figures/1_a.pdf`

**Key Features:**
- Categorizes samples by cancer type and histone modification
- Displays sample sizes (N=448 for H3K4me3, N=303 for H3K27ac)
- Color-coded by assay type

### `ctDNA_levels_by_cancertype.R`

Creates boxplots with jittered points showing ctDNA levels across different cancer types.

**Input Files:**
- Sample metadata Excel file
- Sample information files

**Output:**
- Boxplot with jittered points showing ctDNA distribution
- PDF figure saved to `Manuscripts/Figures/1_b.pdf`

**Key Features:**
- Displays median, quartiles, and individual data points
- Sorted by median ctDNA levels
- Color-coded visualization

### `tissue_overlap_difference_barplot.R`

Generates horizontal bar plots comparing tissue overlap percentages between cfChIP and WBC ChIP data across different Roadmap Epigenomics tissues.

**Input Files:**
- Combined results file from EpiMap overlap analysis
- Chromatin state information (e.g., `EnhA1_EnhA2_TssA_EnhG1_EnhG2`)

**Output:**
- Horizontal bar plot showing difference in overlap percentages
- CSV tables for supplementary data (sup1.csv, sup3.csv)
- PDF figure saved to `Manuscripts/Figures/1_c_*.pdf`

**Key Features:**
- Calculates ratio and delta (difference) between cfChIP and WBC ChIP
- Highlights top 15 and bottom 10 tissues by ratio
- Color-coded by direction of difference (cfChIP > WBC or vice versa)

### `tissue_enrichment_pie_chart.R`

Creates pie charts summarizing tissue enrichment patterns, showing the proportion of tissues where cfChIP or WBC ChIP shows higher overlap.

**Input Files:**
- Combined results file from EpiMap overlap analysis

**Output:**
- Pie chart showing enrichment categories
- PDF figure saved to `Manuscripts/Figures/1_d_*.pdf`

**Key Features:**
- Categorizes tissues as "cfChIP-enriched" or "WBC-enriched"
- Displays percentages for each category
- Handles small slices with appropriate labeling

## Dependencies

- R packages: `ggplot2`, `dplyr`, `readxl`, `tidyr`, `forcats`, `scales`, `extrafont`
- Input data files (see individual script descriptions)

## Usage

Each script is self-contained and can be run independently:

```bash
Rscript sample_composition_by_cancertype.R
Rscript ctDNA_levels_by_cancertype.R
Rscript tissue_overlap_difference_barplot.R
Rscript tissue_enrichment_pie_chart.R
```

**Note:** Update file paths in each script to match your local directory structure. Look for `setwd()` calls and file paths at the beginning of each script.

## Output Format

All scripts generate PDF figures with:
- Consistent color schemes
- High-resolution vector graphics
- Standardized figure dimensions (typically 120x100 mm or similar)

