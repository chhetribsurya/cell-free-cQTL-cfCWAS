# TWAS Data Directory

This directory contains transcriptome-wide association study (TWAS) data and results used in the cfCWAS workflow.

## Contents

- **ONCOARRAY.TWAS.txt**: TWAS results from OncoArray study
- **TCGA_HSQ/**: Directory containing TCGA (The Cancer Genome Atlas) heritability analysis results
  - `TCGA-PRAD.TUMOR.hsq`: Heritability estimates for prostate cancer tumor samples
  - `TCGA-PRAD.TUMOR.pos`: Position file for TCGA prostate cancer analysis
  - `TCGA-PRAD.TUMOR.profile`: Profile file with prediction weights
  - `TCGA-PRAD.TUMOR.profile.err`: Error log file

## Purpose

TWAS data is used to:
- Compare cfCWAS results with transcriptome-wide association findings
- Validate chromatin-trait associations against gene expression associations
- Integrate multi-omics association results

## File Formats

- **.hsq files**: Contain heritability estimates and variance components
- **.pos files**: Contain genomic positions and feature annotations
- **.profile files**: Contain prediction weights and model parameters
- **.txt files**: Tab-delimited association results

## Usage

These files are used for comparative analysis and validation. The workflow may reference these files when comparing cfCWAS results to existing TWAS findings or when performing integrated analyses.

## Note

These are example/reference files from the original CWAS workflow. Users should replace or supplement these with data relevant to their specific analysis.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

