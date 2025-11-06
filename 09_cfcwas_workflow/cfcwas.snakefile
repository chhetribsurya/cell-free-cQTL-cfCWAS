"""
Cell-Free Cistrome-Wide Association Studies (cfCWAS) Snakemake Workflow

This workflow performs cell-free cistrome-wide association studies by:
1. Processing cell-free ChIP-seq data to identify chromatin features
2. Mapping QTLs that influence cell-free chromatin activity
3. Associating predicted chromatin activity with GWAS traits

Adapted from CWAS (https://github.com/scbaca/cwas) for cell-free chromatin data.
"""

configfile: "config.yaml"

# Define sample names from config
samples = config['samples']
gwas_traits = config['gwas_traits']

# Reference files
REFERENCE = config['reference_genome']
LDREF = config['ldref_path']
WASP_PATH = config['wasp_path']
GATK3_PATH = config['gatk3_path']

# Output directories
ANALYSIS_DIR = "analysis"
QTL_DIR = f"{ANALYSIS_DIR}/qtl"
HSQ_DIR = f"{ANALYSIS_DIR}/hsq"
STRATAS_DIR = f"{ANALYSIS_DIR}/stratAS"
FUSION_DIR = f"{ANALYSIS_DIR}/fusion"
CWAS_RESULTS_DIR = f"{ANALYSIS_DIR}/cwas_results"

rule all:
    """
    Main entry point for the workflow.
    """
    input:
        expand(f"{CWAS_RESULTS_DIR}/{{trait}}/cwas.sig.best.txt", trait=gwas_traits),
        expand(f"{QTL_DIR}/{{sample}}/qtl_results.txt", sample=samples)

# ============================================================================
# Quality Control and Preprocessing
# ============================================================================

rule process_cfchip_bam:
    """
    Process cell-free ChIP-seq BAM files using WASP to remove mapping bias.
    """
    input:
        bam = lambda wildcards: config['samples'][wildcards.sample]['bam'],
        vcf = lambda wildcards: config['samples'][wildcards.sample]['vcf']
    output:
        processed_bam = f"{ANALYSIS_DIR}/processed_bam/{{sample}}.wasp.bam",
        stats = f"{ANALYSIS_DIR}/processed_bam/{{sample}}.wasp.stats.txt"
    shell:
        """
        python {WASP_PATH}/mapping/find_intersecting_snps.py \
            --is_paired_end \
            --is_sorted \
            --output_dir {ANALYSIS_DIR}/wasp_temp/{wildcards.sample} \
            --snp_tab {WASP_PATH}/mapping/snp_tab.h5 \
            --snp_index {WASP_PATH}/mapping/snp_index.h5 \
            --haplotype {WASP_PATH}/mapping/haps.h5 \
            --samples {config['samples'][wildcards.sample]['sample_id']} \
            {input.bam} 2> {output.stats}
        """

rule call_peaks:
    """
    Call peaks from cell-free ChIP-seq data.
    """
    input:
        bam = f"{ANALYSIS_DIR}/processed_bam/{{sample}}.wasp.bam"
    output:
        peaks = f"{ANALYSIS_DIR}/peaks/{{sample}}/{{sample}}.peaks.bed"
    params:
        peak_type = lambda wildcards: config['samples'][wildcards.sample].get('peak_type', 'broad')
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -n {wildcards.sample} \
            --outdir {ANALYSIS_DIR}/peaks/{wildcards.sample} \
            -g hs \
            --{params.peak_type} \
            -q 0.05
        """

rule create_consensus_peaks:
    """
    Create consensus peak set from all cell-free ChIP-seq samples.
    """
    input:
        peaks = expand(f"{ANALYSIS_DIR}/peaks/{{sample}}/{{sample}}.peaks.bed", sample=samples)
    output:
        consensus = f"{ANALYSIS_DIR}/peaks/consensus/consensus_peaks.bed"
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | \
        bedtools merge -i - > {output}
        """

# ============================================================================
# QTL Mapping
# ============================================================================

rule map_cfcqtls:
    """
    Map cell-free chromatin QTLs (cfcQTLs) by associating genotypes with chromatin activity.
    """
    input:
        bam = f"{ANALYSIS_DIR}/processed_bam/{{sample}}.wasp.bam",
        peaks = f"{ANALYSIS_DIR}/peaks/consensus/consensus_peaks.bed",
        vcf = lambda wildcards: config['samples'][wildcards.sample]['vcf']
    output:
        qtl_results = f"{QTL_DIR}/{{sample}}/qtl_results.txt",
        qtl_summary = f"{QTL_DIR}/{{sample}}/qtl_summary.txt"
    script:
        "scripts/map_cfcqtls.py"

rule aggregate_qtl_results:
    """
    Aggregate QTL results across all samples.
    """
    input:
        qtl_results = expand(f"{QTL_DIR}/{{sample}}/qtl_results.txt", sample=samples)
    output:
        aggregated = f"{QTL_DIR}/aggregated/all_samples_qtl.txt"
    script:
        "scripts/aggregate_qtl.py"

# ============================================================================
# Model Building and Cross-Validation
# ============================================================================

rule build_prediction_model:
    """
    Build predictive models for chromatin activity using FUSION.
    """
    input:
        qtl_results = f"{QTL_DIR}/aggregated/all_samples_qtl.txt",
        peaks = f"{ANALYSIS_DIR}/peaks/consensus/consensus_peaks.bed"
    output:
        model = f"{FUSION_DIR}/models/{{trait}}/{{trait}}.model.RData",
        weights = f"{FUSION_DIR}/models/{{trait}}/{{trait}}.weights.txt"
    params:
        trait = lambda wildcards: wildcards.trait
    script:
        "scripts/build_prediction_model.R"

rule cross_validate_model:
    """
    Perform cross-validation to assess model performance.
    """
    input:
        model = f"{FUSION_DIR}/models/{{trait}}/{{trait}}.model.RData"
    output:
        cv_results = f"{FUSION_DIR}/cv/{{trait}}/cwas.cv.best.txt"
    script:
        "scripts/cross_validate_model.R"

# ============================================================================
# Association Testing
# ============================================================================

rule run_cwas:
    """
    Run cell-free cistrome-wide association study using FUSION.
    """
    input:
        model = f"{FUSION_DIR}/models/{{trait}}/{{trait}}.model.RData",
        weights = f"{FUSION_DIR}/models/{{trait}}/{{trait}}.weights.txt",
        gwas = f"gwas_data/{{trait}}.sumstats.gz"
    output:
        cwas_results = f"{CWAS_RESULTS_DIR}/{{trait}}/cwas.sig.best.txt",
        cwas_all = f"{CWAS_RESULTS_DIR}/{{trait}}/cwas.all.best.txt"
    params:
        ldref = LDREF,
        trait = lambda wildcards: wildcards.trait
    script:
        "scripts/run_cwas.R"

# ============================================================================
# Heritability Analysis
# ============================================================================

rule calculate_heritability:
    """
    Calculate heritability of cell-free chromatin features.
    """
    input:
        qtl_results = f"{QTL_DIR}/aggregated/all_samples_qtl.txt"
    output:
        hsq_results = f"{HSQ_DIR}/heritability_results.txt"
    script:
        "scripts/calculate_heritability.R"

# ============================================================================
# Stratified Analysis
# ============================================================================

rule stratified_association:
    """
    Perform stratified association analysis by chromatin feature type.
    """
    input:
        cwas_results = expand(f"{CWAS_RESULTS_DIR}/{{trait}}/cwas.sig.best.txt", trait=gwas_traits),
        peaks = f"{ANALYSIS_DIR}/peaks/consensus/consensus_peaks.bed"
    output:
        stratified_results = f"{STRATAS_DIR}/stratified_results.txt"
    script:
        "scripts/stratified_association.R"

