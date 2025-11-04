# Cell-Free Chromatin Quantitative Trait Locus (cQTL) Analysis

This repository contains the computational scripts and analysis pipelines used for identifying and characterizing cell-free chromatin quantitative trait loci (cQTLs) from cell-free ChIP-seq data. The analyses include enrichment testing, developmental specificity assessment, and heritability estimation using stratified linkage disequilibrium score regression (S-LDSC).

## Overview

This repository is organized into eight main analysis modules:

1. **Data Preprocessing** (`01_data_preprocessing/`): Scripts for calculating tissue-specific scores and preparing annotation files
2. **Enrichment Analysis** (`02_enrichment_analysis/`): Scripts for performing cQTL enrichment analyses with permutation testing
3. **LDSC Analysis** (`03_ldsc_analysis/`): Scripts for running stratified LDSC analysis
4. **Sample Metadata Visualization** (`04_sample_metadata/`): Scripts for visualizing sample composition and metadata
5. **Genomic Track Visualization** (`05_genomic_tracks/`): Scripts for creating multi-track genomic visualizations
6. **Allelic Imbalance Visualization** (`06_allelic_imbalance/`): Scripts for visualizing allelic imbalance and Q-value analyses
7. **Enrichment Visualization** (`07_enrichment_visualization/`): Scripts for visualizing enrichment analysis results
8. **Manhattan Plot Visualization** (`08_manhattan_plots/`): Scripts for generating Manhattan plots from CWAS results

## Repository Structure

```
.
├── 01_data_preprocessing/           # Data preprocessing and annotation scripts
├── 02_enrichment_analysis/          # cQTL enrichment analysis scripts
├── 03_ldsc_analysis/                 # Stratified LDSC analysis scripts
├── 04_sample_metadata/               # Sample composition and metadata visualization
├── 05_genomic_tracks/                # Multi-track genomic visualizations
├── 06_allelic_imbalance/             # Allelic imbalance and Q-value visualization
├── 07_enrichment_visualization/      # Enrichment analysis visualization
├── 08_manhattan_plots/               # Manhattan plot generation
└── docs/                             # Documentation and methods
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

## Quick Start

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

- **Sample Metadata** (`04_sample_metadata/`): Sample composition and metadata plots
- **Genomic Tracks** (`05_genomic_tracks/`): Multi-track genomic visualizations
- **Allelic Imbalance** (`06_allelic_imbalance/`): Q-value distributions and allelic fraction analyses
- **Enrichment Visualization** (`07_enrichment_visualization/`): Regulatory and eQTL enrichment plots
- **Manhattan Plots** (`08_manhattan_plots/`): Genome-wide association plots

Each directory contains scripts organized by analysis type rather than figure numbers.

## Analysis Workflow

1. **Data Preprocessing**: Calculate tissue-specific chromatin peak scores by overlapping with Roadmap Epigenomics annotations
2. **Enrichment Analysis**: Test for enrichment of cQTLs in developmentally regulated regulatory regions using Fisher's exact tests and permutation testing
3. **LDSC Analysis**: Estimate heritability enrichment using stratified LDSC
4. **Visualization**: Generate publication-quality figures summarizing the results

## Statistical Methods

### Enrichment Testing

Enrichment analyses use Fisher's exact tests to compare overlap frequencies between cQTL sets and size-matched random backgrounds. Permutation testing (100 permutations) is performed to validate findings and provide empirical p-values.

### Random Background Generation

Random background sets are generated by sampling from the full pool of available cQTLs that passed equivalent significance thresholds. This approach maintains comparable statistical power and preserves genomic distribution characteristics.

### Regulatory Region Classification

Regulatory regions are classified into four mutually exclusive specificity categories:
- Fetal-specific: Active exclusively in fetal tissues
- Stem cell-specific: Active exclusively in stem cell epigenomes
- Developmental-specific: Active in fetal and/or stem cell tissues but absent in adult tissues
- Adult-specific: Active exclusively in adult tissues

## Citation

If you use these scripts in your research, please cite the associated publication.

## License

[Specify license here]

## Contact

For questions or issues, please contact [Your Name/Institution].

## Documentation

Detailed methods are documented in `docs/Methods_section.md`. Each analysis module contains its own README with specific usage instructions.

