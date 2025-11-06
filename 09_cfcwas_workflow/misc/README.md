# Miscellaneous Utilities Directory

This directory contains miscellaneous utility files, reference data, and example files used throughout the cfCWAS workflow.

## Contents

### Reference Files
- **hg19.genome**: Genome size file for hg19 reference build
- **hg19.refGene.tss.tsv**: Transcription start site (TSS) coordinates for hg19
- **glist-hg19**: Gene list file for hg19

### Example Data Files
- **LNCaP_parental_H3K27ac.rep1_sorted_peaks.narrowPeak.bed**: Example H3K27ac peak file
- **LNCaP_parental_AR.rep1_sorted_peaks.narrowPeak.bed**: Example AR (Androgen Receptor) peak file
- **LNCaP_FitHiChIP.interactions_FitHiC_Q0.01_IGV.bed**: Example chromatin interaction file

### Utility Files
- **dummysnps.rds**: R data file with dummy SNP data for testing
- **Nanog.ctrl.motif**: Example motif file
- **predict.head**: Header file for prediction output format

## Purpose

These files serve various purposes:
- Reference genome annotations and coordinates
- Example data for testing and validation
- Format templates for output files
- Utility data structures

## Usage

Most of these files are used internally by workflow scripts. Users typically do not need to modify these files unless customizing the workflow for specific use cases.

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

