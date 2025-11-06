# Cell-Free Cistrome-Wide Association Studies (cfCWAS)

[![Snakemake](https://img.shields.io/badge/Snakemake-6.0%2B-orange.svg)](https://snakemake.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![FUSION](https://img.shields.io/badge/FUSION-1.0%2B-green.svg)](https://github.com/gusevlab/fusion_twas)
[![WASP](https://img.shields.io/badge/WASP-0.3.4%2B-purple.svg)](https://github.com/bmvdgeijn/WASP)

This directory contains the complete workflow for performing cell-free cistrome-wide association studies (cfCWAS). cfCWAS identifies chromatin features from cell-free circulating chromatin that are genetically associated with traits of interest.

## Overview

Cell-free Cistrome-Wide Association Studies (cfCWAS) extend the concept of Cistrome-Wide Association Studies (CWAS) to cell-free chromatin data. This approach:

1. Uses cell-free ChIP-seq (cfChIP-seq) data from multiple individuals to model chromatin activity as a function of nearby SNP genotypes
2. Associates predicted chromatin activity with traits of interest using GWAS summary statistics
3. Identifies genetic variants that influence chromatin accessibility and histone modification patterns in cell-free DNA

## Key Differences from Traditional CWAS

- **Data Source**: Uses cell-free ChIP-seq data (cfChIP-seq) instead of traditional ChIP-seq from cell lines or tissues
- **Sample Characteristics**: Analyzes cell-free circulating chromatin from plasma/serum samples
- **Biological Context**: Captures chromatin states from diverse cell types present in circulation
- **Applications**: Particularly relevant for cancer genomics, liquid biopsies, and circulating biomarkers

## Repository Structure

```
09_cfcwas_workflow/
├── cfcwas.snakefile          # Main Snakemake workflow
├── config.yaml                # Configuration file (update with your data paths)
├── submit.sh                  # SLURM submission script
├── env/
│   └── cfcwas_env.yml         # Conda environment specification
├── scripts/                   # Analysis scripts
├── data/
│   └── peaks/                 # Cell-free chromatin peak files
│       └── consensus/         # Consensus peak sets
├── LDREF/                     # LD reference data from 1000 Genomes
├── gwas_data/                 # GWAS summary statistics files
├── qtl/                       # QTL analysis results
├── hsq/                       # Heritability analysis results
├── stratAS/                   # Stratified association analysis
├── fusion/                    # FUSION software integration
├── twas_data/                 # TWAS-related data
├── run_files/                 # Example configuration files
└── misc/                      # Miscellaneous utilities
```

## Requirements

### Software Dependencies

The workflow requires:

1. **Conda**: For environment management
2. **Snakemake**: For workflow orchestration
3. **Python**: For data processing scripts
4. **R**: For statistical analyses
5. **WASP**: For mitigating mapping bias for reads covering heterozygous SNPs
   - Repository: https://github.com/bmvdgeijn/WASP
6. **GATK3**: For variant calling and processing
   - Download: https://software.broadinstitute.org/gatk/download/archive
7. **FUSION**: For TWAS analysis
   - Reference data: https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2

### Data Requirements

For each cell-free ChIP-seq sample, you need:

1. **BAM files**: Aligned cell-free ChIP-seq reads
2. **BED files**: Peak coordinates from cfChIP-seq data
3. **VCF files**: Phased genotypes (phasing performed with Eagle2 using Sanger Imputation Service: https://imputation.sanger.ac.uk/)

## Installation

### 1. Clone or Download Repository

```bash
cd /path/to/your/repository
```

### 2. Create Conda Environment

```bash
conda env create -f env/cfcwas_env.yml
conda activate cfcwas
```

### 3. Install Additional Software

#### WASP Pipeline

```bash
git clone https://github.com/bmvdgeijn/WASP.git
# Update WASP path in config.yaml
```

#### GATK3

```bash
# Download GATK3 package from Broad Institute
# Extract the .tar file
gatk3-register /path/to/GenomeAnalysisTK.jar
```

#### LDREF Data

```bash
cd LDREF
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar -xjf LDREF.tar.bz2
```

### 4. Configure Workflow

Edit `config.yaml` with your data paths:

```yaml
# Update these paths to match your system
wasp_path: /path/to/WASP
gatk3_path: /path/to/GATK3
ldref_path: /path/to/LDREF
```

### 5. Prepare GWAS Summary Statistics

Place GWAS summary statistics files in `gwas_data/` directory with format:

```
SNP    A1    A2    N      Z
rs123  A     G     100000 2.5
rs456  T     C     100000 -1.8
```

Files should be named `*.sumstats.gz` and gzipped.

## Usage

### Basic Workflow Execution

```bash
# Activate conda environment
conda activate cfcwas

# Run workflow locally (dry run first)
snakemake -n

# Run workflow locally
snakemake --cores 4

# Run on SLURM cluster
sbatch submit.sh
```

### Configuration

The main configuration file is `config.yaml`. Key parameters:

- **Sample information**: BAM files, BED files, VCF files for each sample
- **Reference genome**: Path to reference genome FASTA file
- **Software paths**: WASP, GATK3, FUSION paths
- **Analysis parameters**: Window sizes, significance thresholds, etc.

Example configuration files are provided in `run_files/` directory.

## Workflow Steps

The cfCWAS workflow performs the following steps:

1. **Read Processing**: 
   - Alignment of cell-free ChIP-seq reads
   - Removal of mapping bias using WASP pipeline
   - Quality filtering

2. **Peak Calling**:
   - Identification of chromatin peaks from cfChIP-seq data
   - Consensus peak set generation

3. **QTL Mapping**:
   - Association testing between genotypes and chromatin activity
   - Identification of cell-free chromatin QTLs (cfcQTLs)

4. **Model Building**:
   - Building predictive models for chromatin activity
   - Cross-validation and model selection

5. **Association Testing**:
   - Testing predicted chromatin activity against GWAS traits
   - Identification of significant associations

6. **Heritability Analysis**:
   - Partitioning heritability across chromatin features
   - Estimating contribution of cell-free chromatin to trait heritability

## Output Files

Results are written to the `analysis/` directory:

- `qtl/`: QTL mapping results
- `hsq/`: Heritability analysis results
- `stratAS/`: Stratified association analysis
- `fusion/`: FUSION/TWAS results
- `cwas_results/`: Final cfCWAS association results

## Key Features

### Cell-Free Chromatin Specific Considerations

1. **Low Input Material**: Workflow optimized for low-input cell-free DNA samples
2. **Multiple Cell Types**: Captures chromatin from diverse circulating cell types
3. **Cancer Applications**: Particularly suited for cancer genomics and liquid biopsy applications
4. **Tissue-Specific Signals**: Can identify tissue-of-origin signals in mixed cell populations

### Statistical Methods

- **QTL Mapping**: Linear mixed models for association testing
- **Permutation Testing**: Empirical significance assessment
- **Cross-Validation**: Model performance evaluation
- **Multiple Testing Correction**: FDR and Bonferroni corrections

## Citation

If you use this workflow, please cite:

1. The original CWAS method: [https://www.nature.com/articles/s41588-022-01168-y]
2. This cfCWAS implementation (when published)

## Troubleshooting

### Common Issues

1. **Memory errors**: Increase memory allocation in `submit.sh` or `config.yaml`
2. **Path errors**: Ensure all paths in `config.yaml` are absolute paths
3. **Missing dependencies**: Verify conda environment is properly activated
4. **VCF format issues**: Ensure VCF files are properly phased and indexed

### Getting Help

For issues specific to:
- **Snakemake**: Check Snakemake documentation
- **WASP**: Check WASP GitHub repository
- **FUSION**: Check FUSION documentation
- **cfCWAS workflow**: Open an issue in the repository

## References

- CWAS methodology: https://www.nature.com/articles/s41588-022-01168-y
- WASP pipeline: https://github.com/bmvdgeijn/WASP
- FUSION software: https://github.com/gusevlab/fusion_twas
- Sanger Imputation Service: https://imputation.sanger.ac.uk/

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

