# Example Configuration Files

[![YAML](https://img.shields.io/badge/YAML-1.2-blue.svg)](https://yaml.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-6.0%2B-orange.svg)](https://snakemake.readthedocs.io/)

This directory contains example configuration files for different analysis scenarios.

## Files

### `config_example_cancer.yaml`

Example configuration for cancer-related cfCWAS analysis.

**Features:**
- Configured for cancer GWAS traits
- Optimized parameters for cancer genomics
- Includes common cancer types (breast, prostate, ovarian, lung, colorectal)

### `config_example_non_cancer.yaml`

Example configuration for non-cancer complex trait analysis.

**Features:**
- Configured for non-cancer GWAS traits
- Includes common complex diseases
- Optimized for broader trait associations

### `config_example_minimal.yaml`

Minimal configuration file for testing the workflow.

**Features:**
- Minimal sample set (1-2 samples)
- Single trait for testing
- Quick execution parameters

## Usage

Copy an example configuration file and modify for your analysis:

```bash
cp run_files/config_example_cancer.yaml config.yaml
# Edit config.yaml with your paths and samples
```

## Configuration Parameters

### Required Parameters

- `wasp_path`: Path to WASP pipeline
- `gatk3_path`: Path to GATK3
- `ldref_path`: Path to LD reference data
- `reference_genome`: Path to reference genome FASTA
- `samples`: Dictionary of sample information
- `gwas_traits`: List of GWAS traits to analyze

### Sample Information

For each sample, provide:
- `bam`: Path to BAM file
- `bed`: Path to BED file (peaks)
- `vcf`: Path to VCF file (phased genotypes)
- `sample_id`: Sample identifier in VCF
- `peak_type`: "broad" or "narrow"

### Optional Parameters

- `qtl`: QTL mapping parameters
- `model`: Model building parameters
- `association`: Association testing parameters
- `peak_calling`: Peak calling parameters
- `cluster`: SLURM cluster configuration

## Notes

- All paths should be absolute paths
- Ensure VCF files are properly phased and indexed
- BAM files should be sorted and indexed
- Peak BED files should be in standard BED format (chr, start, end)

## Contact

For questions or issues, please contact: surya_chhetri@dfci.harvard.edu

