# Genomic Track Visualization

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14%2B-green.svg)](https://bioconductor.org/)
[![Gviz](https://img.shields.io/badge/Gviz-1.40%2B-blue.svg)](https://bioconductor.org/packages/Gviz/)
[![GenomicRanges](https://img.shields.io/badge/GenomicRanges-1.46%2B-green.svg)](https://bioconductor.org/packages/GenomicRanges/)

This directory contains R scripts for creating genomic track visualizations using the Gviz package. These scripts generate publication-quality multi-track plots showing chromatin signals, gene annotations, GWAS data, and SNP positions at specific genomic loci.

## Scripts

### `genomic_track_locus_chr5.R`

Creates a multi-track genomic visualization for a specific locus on chromosome 5, showing cfChIP-seq and WBC ChIP-seq signals, GWAS data, gene annotations, and SNP positions.

**Genomic Region:** chr5:1889346 (with ±4kb window)

**Tracks Included:**
- Genome axis track
- GWAS significance track (-log10 p-values)
- cfChIP-seq H3K27ac signal track
- WBC ChIP-seq H3K27ac signal track
- SNP/indel annotation track
- Gene annotation track (with gene symbols)

**Input Files:**
- BigWig files for cfChIP and WBC ChIP signals
- GWAS summary statistics file (PRAD.GWAS.sumstats.txt)
- Reference genome annotations (hg19)

**Output:**
- PDF figure saved to `Manuscripts/Figures/1_f_raw.pdf`

### `genomic_track_locus_chr17_HOXB13.R`

Creates a multi-track genomic visualization for the HOXB13 gene region on chromosome 17, showing signals from multiple cancer types (NEPC, CRC, BRCA, PRAD) and WBC.

**Genomic Region:** chr17:46802126-46806112 (HOXB13 gene)

**Tracks Included:**
- Multiple cancer type-specific cfChIP-seq tracks
- WBC ChIP-seq track
- Gene annotation track
- Genome axis track

**Input Files:**
- BigWig files for different cancer types (NEPC, CRC, BRCA, PRAD)
- WBC BigWig files
- Reference genome annotations (hg19)

**Output:**
- PDF figure showing chromatin signals across cancer types

### `genomic_track_locus_chr17_general.R`

Creates a multi-track genomic visualization for a general locus on chromosome 17, showing cfChIP-seq and WBC ChIP-seq signals with GWAS data.

**Genomic Region:** chr8:128901991-128901992 (with ±15kb window)

**Tracks Included:**
- GWAS significance track
- cfChIP-seq H3K27ac signal track
- WBC ChIP-seq H3K27ac signal track
- Gene annotation track
- Genome axis track

**Input Files:**
- BigWig files for cfChIP and WBC ChIP signals
- GWAS summary statistics file (RCC.GCST90320055.GWAS.chr8.sumstats.txt)
- Reference genome annotations (hg19)

**Output:**
- PDF figure saved to `Manuscripts/Figures/1_g_raw.pdf`

## Dependencies

- R packages: `Gviz`, `GenomicRanges`, `rtracklayer`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `Homo.sapiens`, `data.table`, `dplyr`, `ggplot2`, `extrafont`
- Bioconductor packages: Available through `BiocManager::install()`

## Installation

Install required Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer", 
                       "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                       "Homo.sapiens", "AnnotationDbi"))
```

## Usage

Each script is self-contained and can be run independently:

```bash
Rscript genomic_track_locus_chr5.R
Rscript genomic_track_locus_chr17_HOXB13.R
Rscript genomic_track_locus_chr17_general.R
```

**Note:** 
- Update file paths in each script to match your local directory structure
- Ensure BigWig files are accessible at the specified paths
- The scripts require a helper functions file (`../Scripts/functions.R`) for track creation utilities

## Customization

### Changing Genomic Regions

Modify the `loc` variable at the beginning of each script:
```r
loc <- "chr5:1889346-1889346"  # Format: "chrX:start-end"
```

### Adjusting Window Size

Modify the `win` variable to change the plotting window:
```r
win = 4e3  # 4kb window
```

### Color Schemes

Each script defines color schemes for different tracks. Modify color variables as needed:
```r
gwas_color <- "#D73027"
cfchip_color <- "blue"
h3k27ac_color <- "grey"
```

## Output Format

All scripts generate publication-ready PDF figures with:
- Helvetica font family
- High-resolution vector graphics
- Standardized track heights and spacing
- Multi-track genomic visualization

## Notes

- All coordinates are in hg19/GRCh37
- BigWig files must be properly indexed
- Global maximum values are calculated across all tracks for consistent y-axis scaling
- Gene symbols are automatically mapped using the Homo.sapiens annotation package

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

