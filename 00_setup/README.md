# Setup and Installation

[![Conda](https://img.shields.io/badge/Conda-4.10%2B-green.svg)](https://docs.conda.io/)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)

**GitHub Repository**: [https://github.com/chhetribsurya/cell-free-cQTL-cfCWAS.git](https://github.com/chhetribsurya/cell-free-cQTL-cfCWAS.git)

This directory provides installation instructions and guidance for configuring the computational environment required for cell-free chromatin QTL analysis.

## Overview

The analysis pipeline requires several software dependencies and tools. This directory provides resources for setting up the necessary computational environment.

## System Requirements

### Operating System
- Linux (recommended) or macOS
- Sufficient disk space for data files and intermediate results
- Access to a compute cluster (recommended for large-scale analyses)

### Software Dependencies

#### Required Tools
- **BEDTools** (v2.30.0 or higher): For genomic interval operations
- **LDSC**: For stratified linkage disequilibrium score regression
- **Python** (3.7+): For analysis scripts
- **R** (4.0+): For statistical analysis and visualization
- **Conda/Mamba**: For environment management
- **nullranges** (Bioconductor): For covariate-matched genomic ranges (`matchRanges`)
- **UCSC kent utilities**: `bigWigAverageOverBed` for mappability / accessibility covariates

#### Optional but Recommended
- **SLURM**: For job scheduling on compute clusters
- **WASP**: For mapping bias correction in allelic imbalance analysis
- **FUSION**: For cistrome-wide association studies
- **OpenSSL**: For reproducible `shuf --random-source` streams in bootstrap analyses

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/chhetribsurya/cell-free-cQTL-cfCWAS.git
cd cell-free-cQTL-cfCWAS
```

### 2. Set Up Python Environment

```bash
# Create conda environment for Python dependencies
conda env create -f 09_cfcwas_workflow/env/cfcwas_env.yml
conda activate cfcwas
```

### 3. Install R Packages

```r
# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "data.table", 
                   "RColorBrewer", "readr", "scales", "patchwork",
                   "ggpubr", "extrafont"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "Homo.sapiens", "AnnotationDbi", "qvalue",
                       "ComplexHeatmap", "circlize", "nullranges"))
```

### 4. Install External Tools

#### BEDTools
```bash
# Using conda (recommended)
conda install -c bioconda bedtools

# Or compile from source
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
tar -zxvf bedtools-2.30.0.tar.gz
cd bedtools2
make
```

#### LDSC
```bash
# Clone LDSC repository
git clone https://github.com/bulik/ldsc.git
cd ldsc
# Follow LDSC installation instructions
# Note: LDSC requires Python 2.7
```

#### UCSC kent utilities (`bigWigAverageOverBed`)
```bash
# Using conda (recommended)
conda install -c bioconda ucsc-bigwigaverageoverbed

# Required by 10_matched_random_background/scripts/make_matched_random_bg.sh
```

### 5. Download Reference Data

#### Reference Genome
- Download hg19/GRCh37 reference genome FASTA files
- Index using samtools: `samtools faidx reference.fa`

#### Roadmap Epigenomics Data
- Download chromatin state annotations from Roadmap Epigenomics Project
- Organize by reference epigenome

#### 1000 Genomes Project Data
- Download reference panels for LDSC analysis
- Download HapMap3 SNP lists

## Configuration

### Environment Variables

Set the following environment variables in your shell configuration (e.g., `~/.bashrc` or `~/.zshrc`):

```bash
export LDSC_PATH=/path/to/ldsc
export BEDTOOLS_PATH=/path/to/bedtools
export REFERENCE_GENOME=/path/to/reference/genome
export ROADMAP_DATA=/path/to/roadmap/epigenomics/data
```

### Path Configuration

Update path references in scripts:
- Check `01_data_preprocessing/` scripts for data directory paths
- Update `config.yaml` in `09_cfcwas_workflow/` with your data paths
- Update resource paths in `10_matched_random_background/run_code.sh`
- Update cluster data paths in `11_roadmap_bootstrap_loo/scripts/`
- Modify script-specific paths as needed

## Verification

Test your installation:

```bash
# Test BEDTools
bedtools --version

# Test Python environment
python -c "import pandas, numpy, pybedtools; print('Python packages OK')"

# Test R packages
Rscript -e "library(ggplot2); library(GenomicRanges); cat('R packages OK\n')"
```

## Troubleshooting

### Common Issues

1. **BEDTools not found**: Ensure BEDTools is in your PATH or update script paths
2. **Python package errors**: Verify conda environment is activated
3. **R package installation fails**: Check Bioconductor version compatibility
4. **Permission errors**: Ensure write permissions for output directories

### Getting Help

For installation issues:
- Check individual module READMEs for specific requirements
- Review error messages for missing dependencies
- Verify all paths are correctly configured

## Next Steps

After completing setup:
1. Review the main `README.md` for overview
2. Start with `01_data_preprocessing/` for data preparation
3. Follow the Quick Start guide in the main README

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

