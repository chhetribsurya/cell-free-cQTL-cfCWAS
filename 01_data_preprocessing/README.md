# Data Preprocessing

This directory contains scripts for preprocessing chromatin peak data and calculating tissue-specific scores.

## Scripts

### `calculate_tissue_specific_score.sh`

Calculates tissue-specific scores for chromatin peaks by overlapping them with Roadmap Epigenomics Project chromatin state annotations across multiple reference epigenomes.

**Usage:**
```bash
bash calculate_tissue_specific_score.sh <query_peaks.bed> <assay_name>
```

**Parameters:**
- `query_peaks.bed`: Input BED file containing chromatin peaks (e.g., `data/cfChIP_H3K27ac/consensus.peaks.bed`)
- `assay_name`: Name identifier for the assay (e.g., `cfChIP.H3K27ac`)

**Output:**
- `<assay_name>.peak.overlaped.tsv`: Tab-separated file with overlap information for each peak across all reference epigenomes
- `<assay_name>.peak.freq.bed`: BED file with frequency counts of tissue overlap for each peak

**Requirements:**
- BEDTools (v2.30.0 or higher)
- Roadmap Epigenomics Project chromatin state annotation files
- Reference epigenome BED files in the specified directory structure

**Method:**
The script:
1. Sorts the query peaks by chromosome and start position
2. Intersects each peak with chromatin state annotations from 98 reference epigenomes
3. Generates a binary overlap matrix indicating which peaks overlap with each reference epigenome
4. Calculates the frequency of tissue overlap for each peak

### `plot_tissue_specific_score.R`

R script for visualizing tissue-specific scores calculated from the overlap analysis.

**Usage:**
```bash
Rscript plot_tissue_specific_score.R
```

**Input Files:**
- `Results/EpiMap_overlap/random.peak.freq.bed`: Random background peaks frequency data
- `Results/EpiMap_overlap/cfChIP.K27.peak.freq.bed`: cfChIP H3K27ac peaks frequency data
- `Results/EpiMap_overlap/wbc.peak.freq.bed`: WBC H3K27ac peaks frequency data

**Output:**
- Histogram plots showing the distribution of tissue overlap percentages
- Category breakdowns (tissue-specific, moderately specific, broadly active)

**Method:**
Peaks are categorized based on the percentage of tissues in which they are active:
- No overlap (0%): Not active in any tissue
- Tissue-specific (1-25%): Active in a small subset of tissues
- Moderately specific (25-75%): Active in an intermediate number of tissues
- Broadly active (>75%): Active in most tissues

## Dependencies

- BEDTools v2.30.0
- R (ggplot2, dplyr, data.table, scales, RColorBrewer)

## Notes

- The script assumes a specific directory structure for Roadmap Epigenomics data
- Chromosome coordinates must include the "chr" prefix
- All genomic intervals are in 0-based, half-open BED format

