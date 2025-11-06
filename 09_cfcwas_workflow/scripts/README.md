# cfCWAS Analysis Scripts

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![pybedtools](https://img.shields.io/badge/pybedtools-0.9.0%2B-yellow.svg)](https://daler.github.io/pybedtools/)
[![FUSION](https://img.shields.io/badge/FUSION-1.0%2B-green.svg)](https://github.com/gusevlab/fusion_twas)

This directory contains the core analysis scripts for the cell-free cistrome-wide association studies (cfCWAS) workflow.

## Scripts Overview

### QTL Mapping

#### `map_cfcqtls.py`

Maps cell-free chromatin QTLs (cfcQTLs) by associating SNP genotypes with chromatin activity measured from cell-free ChIP-seq data.

**Purpose:**
- Identifies genetic variants that influence chromatin accessibility and histone modification patterns in cell-free DNA
- Tests associations between nearby SNPs and chromatin peak intensities
- Accounts for cell-free DNA specific characteristics (low input, mixed cell populations)

**Input:**
- Processed BAM files (WASP-processed cell-free ChIP-seq reads)
- Consensus peak BED file
- Phased VCF file with genotypes

**Output:**
- QTL results file with association statistics
- QTL summary statistics

**Key Features:**
- Uses linear mixed models for association testing
- Accounts for read depth and sequencing quality
- Handles low-input cell-free DNA samples
- Corrects for multiple testing

#### `aggregate_qtl.py`

Aggregates QTL results across all cell-free ChIP-seq samples to create a unified QTL map.

**Purpose:**
- Combines QTL findings from multiple samples
- Identifies consistent associations across samples
- Creates meta-analysis QTL results

**Input:**
- Individual sample QTL results files

**Output:**
- Aggregated QTL results file

### Model Building

#### `build_prediction_model.R`

Builds predictive models for chromatin activity using FUSION software, adapted for cell-free chromatin data.

**Purpose:**
- Creates models that predict chromatin activity from genotype data
- Uses cross-validation to select optimal models
- Generates weights for each chromatin feature

**Input:**
- Aggregated QTL results
- Consensus peak set

**Output:**
- Model RData file
- Model weights file

**Key Features:**
- Adapted from FUSION TWAS methodology
- Optimized for cell-free chromatin data characteristics
- Cross-validation for model selection

#### `cross_validate_model.R`

Performs cross-validation to assess model performance and select optimal models.

**Purpose:**
- Evaluates model prediction accuracy
- Selects best models using cross-validation
- Estimates model R-squared and performance metrics

**Input:**
- Model files

**Output:**
- Cross-validation results file (`cwas.cv.best.txt`)

### Association Testing

#### `run_cwas.R`

Performs cell-free cistrome-wide association studies by testing predicted chromatin activity against GWAS traits.

**Purpose:**
- Associates predicted cell-free chromatin activity with complex traits
- Identifies chromatin features genetically associated with disease
- Uses FUSION software for association testing

**Input:**
- Model files
- Model weights
- GWAS summary statistics files

**Output:**
- Significant association results (`cwas.sig.best.txt`)
- All association results (`cwas.all.best.txt`)

**Key Features:**
- Uses FUSION association testing framework
- Accounts for LD structure using reference panels
- Multiple testing correction
- Outputs both significant and all associations

### Heritability Analysis

#### `calculate_heritability.R`

Calculates the heritability of cell-free chromatin features.

**Purpose:**
- Estimates how much of the variation in chromatin activity is due to genetics
- Partitions heritability across different chromatin features
- Quantifies the contribution of cell-free chromatin to trait heritability

**Input:**
- Aggregated QTL results

**Output:**
- Heritability results file

**Key Features:**
- Uses variance component models
- Estimates narrow-sense heritability
- Accounts for cell-free DNA specific characteristics

### Stratified Analysis

#### `stratified_association.R`

Performs stratified association analysis by chromatin feature type, histone modification, or other characteristics.

**Purpose:**
- Identifies which types of chromatin features are most associated with traits
- Compares associations across different histone modifications
- Identifies tissue-specific or cell-type-specific signals

**Input:**
- CWAS association results
- Consensus peak set with annotations

**Output:**
- Stratified association results

**Key Features:**
- Stratifies by peak type, histone mark, chromatin state
- Statistical comparisons between strata
- Visualization-ready output

## Usage

These scripts are typically called by the Snakemake workflow (`cfcwas.snakefile`) and should not be run directly. However, they can be run independently for testing or custom analyses.

### Python Scripts

```bash
# Example: Run QTL mapping for a single sample
python scripts/map_cfcqtls.py \
    --bam analysis/processed_bam/sample1.wasp.bam \
    --peaks analysis/peaks/consensus/consensus_peaks.bed \
    --vcf data/genotypes/sample1.phased.vcf.gz \
    --output analysis/qtl/sample1/
```

### R Scripts

```bash
# Example: Run association testing
Rscript scripts/run_cwas.R \
    --model analysis/fusion/models/trait1/trait1.model.RData \
    --weights analysis/fusion/models/trait1/trait1.weights.txt \
    --gwas gwas_data/trait1.sumstats.gz \
    --output analysis/cwas_results/trait1/
```

## Dependencies

### Python Scripts
- pandas
- numpy
- scipy
- pybedtools
- pysam
- scikit-learn

### R Scripts
- FUSION software (https://github.com/gusevlab/fusion_twas)
- data.table
- dplyr
- Matrix
- Rcpp
- RcppArmadillo

## Notes

- All scripts are designed to work with hg19/GRCh37 reference genome
- VCF files must be phased (use Eagle2 or similar)
- BAM files should be processed through WASP pipeline to remove mapping bias
- Scripts handle low-input cell-free DNA samples appropriately
- Output files follow standard formats for downstream analysis

## Citation

When using these scripts, please cite:
1. The original CWAS methodology
2. FUSION software (for model building and association testing)
3. This cfCWAS implementation (when published)

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

