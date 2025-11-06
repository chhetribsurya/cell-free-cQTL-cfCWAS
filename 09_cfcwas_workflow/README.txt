Cell-Free Cistrome-Wide Association Studies (cfCWAS)

This repository contains scripts and workflows for performing cell-free cistrome-wide association studies (cfCWAS), adapted from the original cistrome-wide association study (CWAS) framework described in https://www.nature.com/articles/s41588-022-01168-y. cfCWAS extends the CWAS methodology to analyze cell-free circulating chromatin from plasma/serum samples, enabling the identification of chromatin features that are genetically associated with traits of interest using liquid biopsy approaches.

OVERVIEW

cfCWAS is a computational framework for identifying genetic variants that influence chromatin accessibility and histone modification patterns in cell-free DNA. Unlike traditional CWAS that analyzes chromatin from cell lines or tissues, cfCWAS leverages cell-free chromatin immunoprecipitation sequencing (cfChIP-seq) data to map genetic effects on regulatory element activity in circulating chromatin. This approach provides unique advantages:

- Access to chromatin states from diverse cell types present in circulation
- Non-invasive sampling through liquid biopsies
- Ability to study chromatin from tissues that are typically inaccessible
- Scalable analysis of large cohorts

The cfCWAS framework identifies chromatin quantitative trait loci (cQTLs) and performs association testing between predicted chromatin activity and disease traits using GWAS summary statistics.

METHODOLOGY

The cfCWAS workflow consists of two main steps:

1. QTL Mapping: Use cell-free ChIP-seq (cfChIP-seq) data from multiple individuals to model chromatin activity as a function of nearby SNP genotypes. This step identifies chromatin QTLs (cQTLs) that influence histone modification patterns (e.g., H3K27ac, H3K4me3) in cell-free chromatin.

2. Association Testing: Associate predicted chromatin activity with traits of interest using GWAS summary statistics. This identifies chromatin features (peaks, regulatory elements) that are genetically associated with disease risk or other complex traits.

The analysis incorporates several key components:
- Allelic imbalance detection to identify variant-specific chromatin activity
- Mapping bias correction using the WASP pipeline
- Prediction weight computation using FUSION
- Stratified association analysis for context-dependent effects

PUBLICATION

This workflow is based on the CWAS methodology described in:
Baca, S.C. et al. Genetic determinants of chromatin reveal prostate cancer risk mediated by context-dependent gene regulation. Nature Genetics 54, 1364-1375 (2022). https://www.nature.com/articles/s41588-022-01168-y

Please see the publication for detailed methodology, validation experiments, and biological interpretation of results.

WORKFLOW IMPLEMENTATION

This analysis is incorporated into a Snakemake workflow in cfcwas.snakefile (adapted from cwas.snakefile). The analysis can be performed in a conda environment whose recipe is provided in env/cfcwas_env.yml.

SETUP AND INSTALLATION

To set up and run the cfCWAS workflow:

1. Environment Setup
   Create and activate the conda environment specified by env/cfcwas_env.yml:
   
   conda env create -f env/cfcwas_env.yml
   conda activate cfcwas

2. Configuration
   Update the config.yaml file with the locations of reference data files. For each sample, you need:
   
   (a) BAM file: Aligned cell-free ChIP-seq reads (cfChIP-seq data)
       - Data should be from H3K27ac or H3K4me3 cfChIP-seq experiments
       - BAM files should be sorted and indexed
   
   (b) BED file: Peak coordinates from the cfChIP-seq data
       - Generated from peak calling on cell-free chromatin data
       - Standard BED format with peak coordinates
   
   (c) VCF file: Phased genotypes for each sample
       - Phased genotypes are required for allelic imbalance analysis
       - For this project, phasing was performed with Eagle2 using the Sanger Imputation Service (https://imputation.sanger.ac.uk/)
       - VCF files should be indexed (tabix)
   
   Example configuration files are provided in run_files/ directory, including:
   - config.H3K27ac.yaml: Example configuration for H3K27ac cfChIP-seq data
   - config.AR.yaml: Example configuration for androgen receptor ChIP-seq data
   - samples.H3K27ac.txt: Sample list format
   - samples.AR.txt: Sample list format

3. Required Data and Software Downloads
   
   (a) WASP Pipeline
       The WASP pipeline is required for mitigating mapping bias for reads covering heterozygous SNPs.
       Repository: https://github.com/bmvdgeijn/WASP
       Clone the repository and specify its location in the config.yaml file.
       This is critical for accurate allelic imbalance detection in cell-free chromatin data.
   
   (b) GATK3 Package
       Download GATK3 package (package-archive_gatk_GenomeAnalysisTK-3.8-0-ge9d806836.tar) from:
       https://software.broadinstitute.org/gatk/download/archive
       Expand the .tar file and activate it using:
       gatk3-register {path to gatk3 jar file}
   
   (c) LD Reference Data
       Download LD reference data from 1000 Genomes Project:
       https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
       Extract the archive into the LDREF/ directory:
       cd LDREF
       wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
       tar -xjf LDREF.tar.bz2
   
   (d) GWAS Summary Statistics
       Place GWAS summary statistics files in the gwas_data/ directory.
       Files should be named with the pattern *.sumstats.gz (gzipped tab-delimited format).
       
       Required format:
       SNP    A1    A2    N      Z
       rs7899632    A    G    140254    -2.98765
       rs3750595    A    C    140254    2.98765
       rs10786405    T    C    140254    2.90123
       
       Where:
       - SNP: SNP identifier (rsID)
       - A1: Effect allele
       - A2: Alternative allele
       - N: Sample size
       - Z: Z-score (effect size / standard error)

4. Activate Environment
   Activate the cfCWAS conda environment:
   conda activate cfcwas

5. Run Workflow
   Execute the Snakemake workflow. For a compute cluster with SLURM workload manager:
   sbatch submit.sh
   
   For local execution:
   snakemake --cores <number_of_cores>
   
   Output files will populate in the analysis/ directory, including:
   - QTL mapping results (qtl/)
   - Heritability estimates (hsq/)
   - Association test results (cwas_results/)
   - Stratified association analysis (stratAS/)

DIRECTORY STRUCTURE

- scripts/: Analysis scripts for QTL mapping, enrichment testing, and visualization
- fusion/: FUSION software integration scripts
- hsq/: Heritability analysis and weight computation scripts
- qtl/: Reference QTL data files (tissue-specific chromatin QTLs)
- stratAS/: Stratified association analysis scripts
- data/: Cell-free chromatin peak data and consensus peak sets
- gwas_data/: GWAS summary statistics files
- LDREF/: Linkage disequilibrium reference data from 1000 Genomes
- run_files/: Example configuration files
- env/: Conda environment specifications
- misc/: Reference files and utilities

KEY DIFFERENCES FROM TRADITIONAL CWAS

- Data Source: Uses cell-free ChIP-seq (cfChIP-seq) instead of traditional ChIP-seq from cell lines
- Sample Type: Analyzes circulating chromatin from plasma/serum samples
- Biological Context: Captures chromatin states from diverse cell types in circulation
- Applications: Particularly relevant for cancer genomics, liquid biopsies, and circulating biomarkers
- Allelic Imbalance: Enhanced focus on detecting allelically imbalanced chromatin activity in cell-free DNA

For detailed usage instructions and additional documentation, see README.md in this directory.
