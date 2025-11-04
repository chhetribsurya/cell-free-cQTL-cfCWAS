#!/usr/bin/env python3
"""
cQTL vs Random cQTL CRE Enrichment with Permutations
====================================================

This script compares input cQTL BED vs random cQTL BED for overlap with a CRE BED file,
performs Fisher's exact test, and runs 100 permutations by sampling the random cQTLs
to match the input count. Outputs summary statistics and permutation results.

Usage:
    python cQTL_CRE_enrichment_permutation.py \
        --input_cqtl input_cqtls_with_id.bed \
        --random_cqtl random_cqtls_with_id.bed \
        --cre_bed ctDNA_CRE_ChIP_peaks_positive.bed \
        --outdir ./enrichment_results

"""
import os
import argparse
import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy import stats
from pathlib import Path
import random
import logging
import matplotlib.pyplot as plt

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="cQTL vs Random cQTL CRE Enrichment with Permutations")
    parser.add_argument('--input_cqtl', required=True, help='Input cQTL BED file')
    parser.add_argument('--random_cqtl', required=True, help='Random cQTL BED file')
    parser.add_argument('--cre_bed', required=True, help='CRE BED file (e.g., ctDNA_CRE_ChIP_peaks_positive.bed)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--n_permutations', type=int, default=100, help='Number of permutations (default: 100)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')
    return parser.parse_args()


def count_overlaps(a_bed, b_bed):
    """Count number of intervals in a_bed that overlap b_bed."""
    return len(a_bed.intersect(b_bed, u=True))


def fisher_test(n_overlap, n_total, n_bg_overlap, n_bg_total):
    table = np.array([
        [n_overlap, n_total - n_overlap],
        [n_bg_overlap, n_bg_total - n_bg_overlap]
    ])
    odds_ratio, p_value = stats.fisher_exact(table)
    return odds_ratio, p_value, table


def odds_ratio_ci(table, alpha=0.05):
    # table: 2x2 numpy array
    # Returns (ci_lower, ci_upper)
    odds_ratio = (table[0,0] * table[1,1]) / (table[0,1] * table[1,0]) if (table[0,1] * table[1,0]) > 0 else float('inf')
    log_or = np.log(odds_ratio) if odds_ratio > 0 else 0
    se_log_or = np.sqrt(np.sum(1 / table))
    z = stats.norm.ppf(1 - alpha/2)
    ci_lower = np.exp(log_or - z * se_log_or)
    ci_upper = np.exp(log_or + z * se_log_or)
    return ci_lower, ci_upper


def process_cqtl_file(cqtl_file, pval_threshold=0.05, qval_threshold=None, genome_build='hg19'):
    """Process cQTL file and convert to BED format.
    
    Args:
        cqtl_file (str): Path to cQTL file
        pval_threshold (float): P-value threshold for filtering (default: 0.05)
        qval_threshold (float): Q-value threshold for filtering (default: None)
        genome_build (str): Genome build to use (default: 'hg19')
    
    Returns:
        pandas.DataFrame: Processed cQTL data in BED format
    """
    logger.info(f"Processing cQTL file: {cqtl_file}")
    
    try:
        # Read cQTL file
        df = pd.read_csv(cqtl_file, sep='\t')
        
        # Filter by p-value or q-value threshold
        original_count = len(df)
        
        # Log the thresholds being used
        logger.info("\nThreshold Information:")
        logger.info(f"  P-value threshold: {pval_threshold}")
        logger.info(f"  Q-value threshold: {qval_threshold}")
        
        if qval_threshold is not None:
            logger.info("\nUsing Q-value based filtering (takes precedence over p-value)")
            logger.info(f"Filtering cQTLs by q-value threshold {qval_threshold}")
            df = df[df['comb.q'] <= qval_threshold]
            filtered_count = len(df)
            logger.info(f"Filtered cQTLs by q-value threshold {qval_threshold}:")
        else:
            logger.info("\nUsing P-value based filtering")
            logger.info(f"Filtering cQTLs by p-value threshold {pval_threshold}")
            df = df[df['pval'] <= pval_threshold]
            filtered_count = len(df)
            logger.info(f"Filtered cQTLs by p-value threshold {pval_threshold}:")
        
        logger.info(f"  Original count: {original_count}")
        logger.info(f"  After filtering: {filtered_count}")
        logger.info(f"  Removed: {original_count - filtered_count} cQTLs")
        
        # Convert to BED format
        df['CHR'] = df['CHR'].astype(str)
        df['CHR'] = df['CHR'].apply(lambda x: x if x.startswith('chr') else 'chr' + x)
        
        # Create BED format
        bed_df = pd.DataFrame({
            'chrom': df['CHR'],
            'start': df['POS'] - 1,  # Convert to 0-based
            'end': df['POS'],
            'identifier': [f"cqtl_{i+1}" for i in range(len(df))]
        })
        
        return bed_df
        
    except Exception as e:
        logger.error(f"Error processing cQTL file: {e}")
        raise


def plot_permutation_histograms(fold_enrichments, odds_ratios, outdir):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.hist(fold_enrichments, bins=20, color='#3498db', alpha=0.7)
    plt.axvline(np.mean(fold_enrichments), color='red', linestyle='--', label=f'Mean: {np.mean(fold_enrichments):.2f}')
    plt.title('Permutation Fold Enrichment')
    plt.xlabel('Fold Enrichment')
    plt.ylabel('Count')
    plt.legend()

    plt.subplot(1, 2, 2)
    finite_odds = [x for x in odds_ratios if np.isfinite(x)]
    if finite_odds:
        plt.hist(finite_odds, bins=20, color='#3498db', alpha=0.7)
        plt.axvline(np.mean(finite_odds), color='red', linestyle='--', label=f'Mean: {np.mean(finite_odds):.2f}')
        plt.title('Permutation Odds Ratio')
        plt.xlabel('Odds Ratio')
        plt.ylabel('Count')
        plt.legend()
    else:
        plt.text(0.5, 0.5, 'No finite odds ratios', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Permutation Odds Ratio')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'permutation_histograms.pdf'), bbox_inches='tight')
    plt.close()


def plot_oddsratio_foldenrich_barplot(odds_ratio, ci_lower, ci_upper, fold_enrichment, outdir):
    fig, ax = plt.subplots(figsize=(6, 6))
    values = [odds_ratio, fold_enrichment]
    bars = ax.bar([0, 1], values,
                  color=['lightcoral', 'orange'], width=0.5, edgecolor='black',
                  label=['Odds Ratio', 'Fold Enrichment'])
    # Add error bar for odds ratio
    ax.errorbar(0, odds_ratio, yerr=[[odds_ratio - ci_lower], [ci_upper - odds_ratio]],
                fmt='none', ecolor='black', capsize=8, lw=2, zorder=10)
    # X labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Odds Ratio', 'Fold Enrichment'])
    # Legends
    from matplotlib.patches import Rectangle
    legend_elements = [Rectangle((0, 0), 1, 1, facecolor='lightcoral', edgecolor='black', label='Odds Ratio'),
                      Rectangle((0, 0), 1, 1, facecolor='orange', edgecolor='black', label='Fold Enrichment')]
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        frameon=False,
        handleheight=1.5,
        handlelength=1.5,
        borderpad=0.8,
        labelspacing=0.8
    )
    # Despine and offset
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_position(('outward', 8))
    ax.spines['bottom'].set_position(('outward', 8))
    ax.tick_params(axis='both', which='both', length=6, width=1, direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_ylabel('Value')
    # Ensure y-axis starts and ends with ticks, rounded to integer if possible
    y_min = 0
    y_max = max(ci_upper, fold_enrichment) * 1.15
    # Choose integer ticks if range is reasonable, else use 2 decimals
    if y_max <= 20:
        yticks = np.arange(np.floor(y_min), np.ceil(y_max)+1, 1)
    else:
        yticks = np.linspace(y_min, y_max, num=6)
        yticks = np.round(yticks).astype(int)
    if yticks[-1] < y_max:
        yticks = np.append(yticks, int(np.ceil(y_max)))
    ax.set_ylim([y_min, y_max])
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(int(y)) for y in yticks])
    # Ensure x-axis ticks at both bars, and round if needed
    ax.set_xticks([0, 1])
    # Add value labels
    # Fold enrichment label
    ax.text(1, fold_enrichment + y_max*0.02, f"{fold_enrichment:.2f}", ha='center', va='bottom', fontsize=12, fontweight='bold')
    # Odds ratio label (on top of error bar)
    odds_label_y = ci_upper + y_max*0.02
    ax.text(0, odds_label_y, f"{odds_ratio:.2f}", ha='center', va='bottom', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'oddsratio_foldenrich_barplot.pdf'), bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Load BED files
    input_bed = BedTool(args.input_cqtl)
    random_bed = BedTool(args.random_cqtl)
    cre_bed = BedTool(args.cre_bed)

    # Save the original input cQTL BED file (all columns)
    input_bed_df = pd.read_csv(args.input_cqtl, sep='\t', header=None)
    input_bed_df.to_csv(os.path.join(args.outdir, 'input_cqtl_original.bed'), sep='\t', header=False, index=False)

    # Find and save overlapping intervals (with all columns from original)
    overlap_bed = input_bed.intersect(cre_bed, u=True)
    overlap_lines = set([str(x).strip() for x in overlap_bed])
    overlap_df = input_bed_df[input_bed_df.apply(lambda row: '\t'.join(map(str, row)), axis=1).isin(overlap_lines)]
    overlap_df.to_csv(os.path.join(args.outdir, 'input_cqtl_overlapping_with_CRE.bed'), sep='\t', header=False, index=False)

    n_input = len(input_bed)
    n_random = len(random_bed)

    # Overlap counts
    n_input_overlap = count_overlaps(input_bed, cre_bed)
    n_random_overlap = count_overlaps(random_bed, cre_bed)

    # Fisher's exact test
    odds_ratio, p_value, table = fisher_test(n_input_overlap, n_input, n_random_overlap, n_random)
    ci_lower, ci_upper = odds_ratio_ci(table)
    fold_enrichment = (n_input_overlap / n_input) / (n_random_overlap / n_random) if n_random_overlap > 0 and n_random > 0 else float('inf')

    print(f"Input cQTLs: {n_input}, Overlapping CREs: {n_input_overlap}")
    print(f"Random cQTLs: {n_random}, Overlapping CREs: {n_random_overlap}")
    print(f"Fisher's exact test: Odds ratio = {odds_ratio:.3f}, p-value = {p_value:.3e}")

    # Permutations
    perm_odds = []
    perm_pvals = []
    perm_overlaps = []
    perm_tables = []
    perm_fold_enrichments = []
    perm_ci_lowers = []
    perm_ci_uppers = []
    perm_input_overlaps = []
    perm_input_totals = []
    perm_random_overlaps = []
    perm_random_totals = []
    random_intervals = list(random_bed)
    for i in range(args.n_permutations):
        sampled = random.sample(random_intervals, n_input)
        sampled_bed = BedTool(sampled)
        n_sampled_overlap = count_overlaps(sampled_bed, cre_bed)
        or_perm, pv_perm, tbl_perm = fisher_test(n_input_overlap, n_input, n_sampled_overlap, n_input)
        ci_lower, ci_upper = odds_ratio_ci(tbl_perm)
        fold_enrichment = (n_input_overlap / n_input) / (n_sampled_overlap / n_input) if n_sampled_overlap > 0 else float('inf')
        perm_odds.append(or_perm)
        perm_pvals.append(pv_perm)
        perm_overlaps.append(n_sampled_overlap)
        perm_tables.append(tbl_perm)
        perm_fold_enrichments.append(fold_enrichment)
        perm_ci_lowers.append(ci_lower)
        perm_ci_uppers.append(ci_upper)
        perm_input_overlaps.append(n_input_overlap)
        perm_input_totals.append(n_input)
        perm_random_overlaps.append(n_sampled_overlap)
        perm_random_totals.append(n_input)

    # Empirical p-value: how many permutations have odds >= observed
    emp_p = (sum([x >= odds_ratio for x in perm_odds]) + 1) / (args.n_permutations + 1)

    # Output summary
    summary = {
        'n_input': n_input,
        'n_input_overlap': n_input_overlap,
        'n_random': n_random,
        'n_random_overlap': n_random_overlap,
        'fold_enrichment': fold_enrichment,
        'odds_ratio': odds_ratio,
        'odds_ratio_ci_lower': ci_lower,
        'odds_ratio_ci_upper': ci_upper,
        'fisher_p_value': p_value,
        'empirical_p_value': emp_p,
        'mean_perm_odds': np.mean(perm_odds),
        'std_perm_odds': np.std(perm_odds),
        'mean_perm_overlaps': np.mean(perm_overlaps),
        'std_perm_overlaps': np.std(perm_overlaps),
        'mean_perm_fold_enrichment': np.mean(perm_fold_enrichments),
        'std_perm_fold_enrichment': np.std(perm_fold_enrichments)
    }
    pd.DataFrame([summary]).to_csv(os.path.join(args.outdir, 'enrichment_summary.csv'), index=False)
    pd.DataFrame({
        'input_overlap': perm_input_overlaps,
        'input_total': perm_input_totals,
        'random_overlap': perm_random_overlaps,
        'random_total': perm_random_totals,
        'fold_enrichment': perm_fold_enrichments,
        'odds_ratio': perm_odds,
        'odds_ratio_ci_lower': perm_ci_lowers,
        'odds_ratio_ci_upper': perm_ci_uppers,
        'fisher_p_value': perm_pvals
    }).to_csv(os.path.join(args.outdir, 'permutation_results.csv'), index=False)

    # Plot permutation histograms
    plot_permutation_histograms(perm_fold_enrichments, perm_odds, args.outdir)

    # After summary is computed and saved
    plot_oddsratio_foldenrich_barplot(odds_ratio, ci_lower, ci_upper, fold_enrichment, args.outdir)

    print(f"\nEmpirical p-value (permuted odds >= observed): {emp_p:.4f}")
    print(f"Results saved to {args.outdir}")


if __name__ == "__main__":
    import sys
    # If '--sanity_check' is passed, run with hardcoded files
    if '--sanity_check' in sys.argv:
        # Remove the flag so argparse doesn't complain
        sys.argv.remove('--sanity_check')
        # Hardcoded test paths (update as needed)
        class Args:
            input_cqtl = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cQTL_H3K4me3_qvalBased_V4/original_cqtls_with_id.bed"
            random_cqtl = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cQTL_H3K4me3_qvalBased_V4/random_cqtls_with_id.bed"
            cre_bed = "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/ctDNA_CRE_ChIP_peaks_positive.bed"
            outdir = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cfcQTLs-CREs_enrichment_V2"
            n_permutations = 100
            seed = 42
        args = Args()
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        random.seed(args.seed)
        np.random.seed(args.seed)
        # Load BED files
        input_bed = BedTool(args.input_cqtl)
        random_bed = BedTool(args.random_cqtl)
        cre_bed = BedTool(args.cre_bed)
        
        # Save the original input cQTL BED file (all columns)
        input_bed_df = pd.read_csv(args.input_cqtl, sep='\t', header=None)
        input_bed_df.to_csv(os.path.join(args.outdir, 'input_cqtl_original.bed'), sep='\t', header=False, index=False)

        # Find and save overlapping intervals (with all columns from original)
        overlap_bed = input_bed.intersect(cre_bed, u=True)
        overlap_lines = set([str(x).strip() for x in overlap_bed])
        overlap_df = input_bed_df[input_bed_df.apply(lambda row: '\t'.join(map(str, row)), axis=1).isin(overlap_lines)]
        overlap_df.to_csv(os.path.join(args.outdir, 'input_cqtl_overlapping_with_CRE.bed'), sep='\t', header=False, index=False)

        n_input = len(input_bed)
        n_random = len(random_bed)
        n_input_overlap = count_overlaps(input_bed, cre_bed)
        n_random_overlap = count_overlaps(random_bed, cre_bed)
        odds_ratio, p_value, table = fisher_test(n_input_overlap, n_input, n_random_overlap, n_random)
        ci_lower, ci_upper = odds_ratio_ci(table)
        fold_enrichment = (n_input_overlap / n_input) / (n_random_overlap / n_random) if n_random_overlap > 0 and n_random > 0 else float('inf')
        print(f"Input cQTLs: {n_input}, Overlapping CREs: {n_input_overlap}")
        print(f"Random cQTLs: {n_random}, Overlapping CREs: {n_random_overlap}")
        print(f"Fisher's exact test: Odds ratio = {odds_ratio:.3f}, p-value = {p_value:.3e}")
        perm_odds = []
        perm_pvals = []
        perm_overlaps = []
        perm_tables = []
        perm_fold_enrichments = []
        perm_ci_lowers = []
        perm_ci_uppers = []
        perm_input_overlaps = []
        perm_input_totals = []
        perm_random_overlaps = []
        perm_random_totals = []
        random_intervals = list(random_bed)
        for i in range(args.n_permutations):
            sampled = random.sample(random_intervals, n_input)
            sampled_bed = BedTool(sampled)
            n_sampled_overlap = count_overlaps(sampled_bed, cre_bed)
            or_perm, pv_perm, tbl_perm = fisher_test(n_input_overlap, n_input, n_sampled_overlap, n_input)
            ci_lower, ci_upper = odds_ratio_ci(tbl_perm)
            fold_enrichment = (n_input_overlap / n_input) / (n_sampled_overlap / n_input) if n_sampled_overlap > 0 else float('inf')
            perm_odds.append(or_perm)
            perm_pvals.append(pv_perm)
            perm_overlaps.append(n_sampled_overlap)
            perm_tables.append(tbl_perm)
            perm_fold_enrichments.append(fold_enrichment)
            perm_ci_lowers.append(ci_lower)
            perm_ci_uppers.append(ci_upper)
            perm_input_overlaps.append(n_input_overlap)
            perm_input_totals.append(n_input)
            perm_random_overlaps.append(n_sampled_overlap)
            perm_random_totals.append(n_input)
        emp_p = (sum([x >= odds_ratio for x in perm_odds]) + 1) / (args.n_permutations + 1)
        summary = {
            'n_input': n_input,
            'n_input_overlap': n_input_overlap,
            'n_random': n_random,
            'n_random_overlap': n_random_overlap,
            'fold_enrichment': fold_enrichment,
            'odds_ratio': odds_ratio,
            'odds_ratio_ci_lower': ci_lower,
            'odds_ratio_ci_upper': ci_upper,
            'fisher_p_value': p_value,
            'empirical_p_value': emp_p,
            'mean_perm_odds': np.mean(perm_odds),
            'std_perm_odds': np.std(perm_odds),
            'mean_perm_overlaps': np.mean(perm_overlaps),
            'std_perm_overlaps': np.std(perm_overlaps),
            'mean_perm_fold_enrichment': np.mean(perm_fold_enrichments),
            'std_perm_fold_enrichment': np.std(perm_fold_enrichments)
        }
        pd.DataFrame([summary]).to_csv(os.path.join(args.outdir, 'enrichment_summary.csv'), index=False)
        pd.DataFrame({
            'input_overlap': perm_input_overlaps,
            'input_total': perm_input_totals,
            'random_overlap': perm_random_overlaps,
            'random_total': perm_random_totals,
            'fold_enrichment': perm_fold_enrichments,
            'odds_ratio': perm_odds,
            'odds_ratio_ci_lower': perm_ci_lowers,
            'odds_ratio_ci_upper': perm_ci_uppers,
            'fisher_p_value': perm_pvals
        }).to_csv(os.path.join(args.outdir, 'permutation_results.csv'), index=False)
        # Plot permutation histograms
        plot_permutation_histograms(perm_fold_enrichments, perm_odds, args.outdir)
        # After summary is computed and saved
        plot_oddsratio_foldenrich_barplot(odds_ratio, ci_lower, ci_upper, fold_enrichment, args.outdir)

        print(f"\nEmpirical p-value (permuted odds >= observed): {emp_p:.4f}")
        print(f"Results saved to {args.outdir}")
    else:
        main()
    
    