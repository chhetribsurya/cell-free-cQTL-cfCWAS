# Data Preprocessing

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![BEDTools](https://img.shields.io/badge/BEDTools-2.30%2B-green.svg)](https://bedtools.readthedocs.io/)
[![Bash](https://img.shields.io/badge/Bash-4.0%2B-lightgrey.svg)](https://www.gnu.org/software/bash/)

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

### `epimap_overlap_analysis.sh`

Calculates overlap between cell-free chromatin peaks and individual EpiMap reference epigenomes. This script processes chromatin state patterns and compares cfChIP and WBC ChIP peak overlaps with Roadmap Epigenomics Project reference epigenomes.

**Usage:**
```bash
bash epimap_overlap_analysis.sh "<chromatin_state_patterns>"
```

**Parameters:**
- `chromatin_state_patterns`: Comma-separated list of chromatin state patterns (e.g., "EnhA1,EnhA2,TssA,EnhG1,EnhG2")

**Input Files:**
- EpiMap epigenome data files (compressed BED files)
- Consensus peaks from cfChIP-seq data
- Consensus peaks from WBC ChIP-seq data
- EpiMap label file with epigenome information

**Output:**
- Overlap BED files for each reference epigenome
- Combined results file (`{patterns}_combined_results.tsv`) with:
  - EID (Epigenome ID)
  - WB_count (Whole blood overlap count)
  - CF_count (Cell-free overlap count)
  - Total_count (Total peaks in reference epigenome)
  - NAME (Epigenome name)

**Method:**
1. Selects and merges chromatin states based on specified patterns
2. Calculates overlap between cfChIP peaks and each reference epigenome
3. Downsamples whole-blood peaks by chromosome to match cfChIP distribution
4. Calculates overlap between downsampled WBC peaks and reference epigenomes
5. Combines results into a single summary file

**Key Features:**
- Handles multiple chromatin state patterns
- Performs chromosome-matched downsampling for fair comparison
- Generates comprehensive overlap statistics

### `grouped_epimap_overlap_analysis.sh`

Calculates overlap between cell-free chromatin peaks and grouped EpiMap categories. Groups reference epigenomes by tissue type (e.g., blood vs non-blood) and calculates overlap statistics.

**Usage:**
```bash
bash grouped_epimap_overlap_analysis.sh
```

**Input Files:**
- EpiMap information CSV file with grouping information
- Individual epigenome BED files
- Consensus peaks from cfChIP-seq data
- Consensus peaks from WBC ChIP-seq data

**Output:**
- Grouped merged BED files for each tissue category
- Grouped output ratios file (`EnhA1_EnhA2_EnhG1_EnhG2_grouped_output_ratios.txt`) with:
  - Group name
  - Overlap ratio (percentage)
  - Assay type (cfChIP_H3K27ac or whole-blood)

**Method:**
1. Extracts unique groups from EpiMap information file
2. Concatenates and merges peaks for each group
3. Calculates overlap ratios for cfChIP and WBC data
4. Generates summary file with grouped statistics

**Key Features:**
- Groups epigenomes by tissue categories
- Merges peaks within groups with threshold
- Calculates overlap percentages for each group
- Compares cfChIP and WBC ChIP overlap patterns

### `epimap_overlap_group_comparison.R`

R script for visualizing EpiMap overlap comparisons between blood and non-blood tissue categories. This script generates Extended Figure 1c by comparing overlap frequencies between cfChIP and WBC ChIP data.

**Usage:**
```bash
Rscript epimap_overlap_group_comparison.R
```

**Input Files:**
- EpiMap overlap TSV files (from `epimap_overlap_analysis.sh`):
  - `Results/EpiMap_overlap/wbc.peak.overlaped.tsv`
  - `Results/EpiMap_overlap/cfChIP_H3K27ac.peak.overlaped.tsv`
- EpiMap information CSV file with grouping information

**Output:**
- Bar plot showing mean number of chromatin peaks overlapped by blood vs non-blood categories
- PDF figure saved to `Figures/EpiMap Overlap/NB_B_group_overlapped_mean_{threshold}.pdf`

**Key Features:**
- Calculates overlap counts for top N% of peaks (default: top 25%)
- Groups epigenomes by blood vs non-blood categories
- Performs statistical comparisons (Wilcoxon test) between categories
- Displays mean differences and p-values
- Color-coded by assay type (cfChIP vs WBC ChIP)

**Method:**
1. Reads overlap data from both cfChIP and WBC ChIP
2. Filters to top N% of peaks with highest overlap frequency
3. Counts overlaps per epigenome (EID)
4. Groups by blood/non-blood categories
5. Calculates mean overlap counts and statistical tests
6. Generates bar plot with statistical annotations

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

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

