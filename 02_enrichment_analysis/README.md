# Enrichment Analysis

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![BEDTools](https://img.shields.io/badge/BEDTools-2.30%2B-green.svg)](https://bedtools.readthedocs.io/)
[![pybedtools](https://img.shields.io/badge/pybedtools-0.9.0%2B-yellow.svg)](https://daler.github.io/pybedtools/)

This directory contains scripts for performing enrichment analyses of cell-free chromatin quantitative trait loci (cQTLs) in developmentally regulated regulatory regions and other genomic annotations.

## Scripts

### `cqtl_enrichment_analysis.py`

Main script for analyzing enrichment of cQTLs in developmentally regulated regulatory regions. Performs Fisher's exact tests and permutation testing to assess statistical significance.

**Usage:**
```bash
python cqtl_enrichment_analysis.py \
    --cqtl_file <cqtl_file.bed> \
    --random_background <random_background.bed> \
    --regulatory_regions <regulatory_regions.bed> \
    --output_dir <output_directory> \
    [--n_permutations 100] \
    [--seed 50] \
    [--genome_build hg19]
```

**Parameters:**
- `--cqtl_file`: Input BED file containing cQTLs
- `--random_background`: BED file containing random background cQTLs for comparison
- `--regulatory_regions`: BED file containing regulatory region annotations
- `--output_dir`: Directory for output files
- `--n_permutations`: Number of permutations for empirical p-value calculation (default: 100)
- `--seed`: Random seed for reproducibility (default: 50)
- `--genome_build`: Genome build (hg19 or hg38, default: hg19)

**Output:**
- Enrichment statistics including fold enrichment, odds ratios, and p-values
- Permutation results for empirical significance assessment
- Summary plots and tables

**Method:**
1. Overlaps cQTLs with regulatory region annotations using BEDTools
2. Calculates overlap frequencies for test and background sets
3. Performs Fisher's exact test to assess enrichment
4. Runs permutation testing to validate findings
5. Computes fold enrichment and 95% confidence intervals for odds ratios

### `cfpeak_enrichment_analysis.py`

Analyzes enrichment of cell-free chromatin peaks in regulatory regions. Similar to cQTL enrichment analysis but focused on chromatin peak data.

**Usage:**
```bash
python cfpeak_enrichment_analysis.py \
    --cfpeak_file <cfpeak_file.bed> \
    --random_background <random_background.bed> \
    --regulatory_regions <regulatory_regions.bed> \
    --output_dir <output_directory>
```

**Parameters:** Same as `cqtl_enrichment_analysis.py`

### `cwas_regulatory_intersection.py`

Analyzes the intersection between cfCWAS (cell-free Cistrome-Wide Association Study) loci and 
developmentally regulated regulatory elements. This script performs comprehensive intersection 
analysis, pleiotropy assessment, and trait association analysis for cfCWAS results.

**Usage:**
```bash
python cwas_regulatory_intersection.py
```

**Configuration:**
Update the `analysis_configs` list in the `main()` function with your data paths:
```python
analysis_configs = [
    {
        'name': 'cQTL_cfCWAS_combined',
        'regulatory_dir': "/path/to/regulatory/elements/directory",
        'cwas_file': "/path/to/cfCWAS/results.csv",
        'output_dir': "/path/to/output/directory"
    }
]
```

**Input Requirements:**
- `cwas_file`: CSV file with cfCWAS results containing columns:
  - `ID`: Locus identifier in format "chr:start-end"
  - `pleiotropy`: Pleiotropy classification ("trait-specific" or "pleiotropic")
  - `all_descriptions`: Semicolon-separated trait descriptions
- `regulatory_dir`: Directory containing regulatory element BED files:
  - `adult_specific_intersected.bed`
  - `developmental_specific_intersected.bed`
  - `fetal_specific_intersected.bed`
  - `stem_cell_specific_intersected.bed`

**Output:**
- Intersection statistics and detailed results
- Pleiotropy distribution plots
- Trait association analysis plots (individual and combined)
- Detailed trait analysis files (CSV and TSV)
- Summary report

**Method:**
1. Loads cfCWAS loci and converts to BED format
2. Intersects with regulatory element annotations using pybedtools
3. Analyzes pleiotropy distribution for each regulatory type
4. Performs trait association analysis
5. Generates comprehensive visualizations and reports

### `enrichment_statistics_plots.R`

R script for generating publication-ready plots from enrichment analysis statistics. Creates bar plots for fold enrichment, odds ratios, and fraction overlap across regulatory element categories.

**Usage:**
```r
# Update base_dir and regulatory_element variables at the top of the script
# Then run:
Rscript enrichment_statistics_plots.R
```

**Configuration:**
Update the following variables at the top of the script:
- `regulatory_element`: One of "active", "enhancer", "promoter", or "bivalent"
- `base_dir`: Path to the directory containing enrichment analysis results

**Input Requirements:**
- `detailed_statistics_full_random.csv`: File containing enrichment statistics with columns:
  - `category`: Regulatory category (fetal specific, developmental specific, stem cell specific, adult specific)
  - `fold_enrichment`: Fold enrichment values
  - `odds_ratio`: Odds ratio values
  - `odds_ratio_ci_lower`, `odds_ratio_ci_upper`: Confidence interval bounds
  - `percentage`: Percentage overlap
  - `overlapping_peaks`: Number of overlapping peaks

**Output:**
- `fold_enrichment_plot.pdf`: Bar plot of fold enrichment by category
- `odds_ratio_plot.pdf`: Bar plot of odds ratios with confidence intervals
- `fraction_overlap_plot.pdf`: Bar plot of percentage overlap by category

**Key Features:**
- Color-coded bars by regulatory category:
  - Fetal specific: Yellow (#F1C40F)
  - Developmental specific: Red (#E74C3C)
  - Stem cell specific: Blue (#3498DB)
  - Adult specific: Green (#2ECC71)
- Sorted data by metric for clear visualization
- Error bars for odds ratios showing confidence intervals
- Publication-ready formatting with proper margins and labels

**Dependencies:**
- R (>= 3.6.0)
- ggplot2
- RColorBrewer
- gridExtra

### `cre_enrichment_analysis.py`

Analyzes enrichment of cQTLs within experimentally validated cis-regulatory elements (CREs).

**Usage:**
```bash
python cre_enrichment_analysis.py \
    --input_cqtl <cqtl_file.bed> \
    --random_cqtl <random_cqtl_file.bed> \
    --cre_bed <cre_annotations.bed> \
    --outdir <output_directory> \
    [--n_permutations 100] \
    [--seed 42]
```

**Parameters:**
- `--input_cqtl`: Input cQTL BED file
- `--random_cqtl`: Random cQTL BED file for comparison
- `--cre_bed`: CRE annotation BED file (e.g., `ctDNA_CRE_ChIP_peaks_positive.bed`)
- `--outdir`: Output directory
- `--n_permutations`: Number of permutations (default: 100)
- `--seed`: Random seed (default: 42)

**Output:**
- Summary statistics including overlap counts, fold enrichment, odds ratios, and p-values
- Permutation results

### `qtl_enrichment_analysis.sh`

Shell script wrapper for running enrichment analyses in a batch processing environment (SLURM).

**Usage:**
```bash
bash qtl_enrichment_analysis.sh <foreground_file> <background_file> <line_name> <target_file> <output_directory>
```

**Parameters:**
- `foreground_file`: Foreground cQTL file
- `background_file`: Background cQTL file
- `line_name`: Identifier for the analysis
- `target_file`: Target annotation file
- `output_directory`: Output directory

**Method:**
1. Calculates foreground enrichment
2. Runs 50 permutations on background set
3. Generates random shuffles for additional validation
4. Summarizes results in `enrichment.stats.txt`

**Output:**
- `fg.enrich.txt`: Foreground enrichment results
- Permutation results in `perm/` subdirectory
- `enrichment.stats.txt`: Summary statistics

## Statistical Methods

### Fisher's Exact Test

Enrichment is assessed using Fisher's exact test with a 2×2 contingency table:
- Rows: Test set vs. background set
- Columns: Overlapping vs. non-overlapping regions

### Odds Ratio Calculation

Odds ratios are calculated directly from the contingency table. 95% confidence intervals are computed using the log transformation method (Wald method).

### Permutation Testing

Permutation testing is performed to validate enrichment findings:
- 100 independent random samples are generated from the background set
- Each sample is matched in size to the test set
- Empirical p-values are calculated as: `(n_permutations ≥ observed + 1) / (n_permutations + 1)`

### Fold Enrichment

Fold enrichment is computed as:
```
fold_enrichment = (n_overlap_test / n_total_test) / (n_overlap_random / n_total_random)
```

## Regulatory Region Categories

Regulatory regions are classified into four mutually exclusive categories:
1. **Fetal-specific**: Active exclusively in fetal tissues (13 epigenomes)
2. **Stem cell-specific**: Active exclusively in stem cell epigenomes (22 epigenomes)
3. **Developmental-specific**: Active in fetal and/or stem cell tissues but absent in adult tissues
4. **Adult-specific**: Active exclusively in adult tissues (80 epigenomes)

## Dependencies

- Python 3.x
- pandas
- numpy
- scipy
- pybedtools
- matplotlib
- seaborn
- tqdm
- BEDTools v2.30.0

## Notes

- All genomic coordinates must be in BED format (0-based, half-open)
- Chromosome coordinates must include the "chr" prefix
- Random background sets should be size-matched to test sets for valid comparisons
- Fixed random seeds are used for reproducibility

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

