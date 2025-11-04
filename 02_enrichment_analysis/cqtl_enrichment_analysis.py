#!/usr/bin/env python3
"""
cQTL Analysis with Custom Random Background
==========================================

Analysis pipeline for identifying cell-free circulating chromatin quantitative trait loci (cQTLs)
associated with developmentally regulated regulatory regions using custom random background.

Author: Surya B. Chhetri
"""

import os
import logging
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import random
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import pybedtools
from pybedtools import BedTool
import tempfile
import json
import gzip
import argparse
from collections import OrderedDict
from itertools import combinations
from upsetplot import UpSet, from_memberships
from matplotlib_venn import venn2, venn3
import matplotlib as mpl
from matplotlib import cm


# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress pybedtools warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pybedtools")

def get_chromosome_sizes(genome_build='hg19'):
    """Get chromosome sizes for the specified genome build."""
    if genome_build == 'hg19':
        return {
            'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
            'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
            'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431,
            'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
            'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
            'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
            'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895,
            'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566
        }
    elif genome_build == 'hg38':
        return {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
            'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
            'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
            'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
            'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
            'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
            'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }
    else:
        raise ValueError(f"Unsupported genome build: {genome_build}")


def get_regulatory_states(state_type='enhancer'):
    """Get regulatory state definitions based on the specified type.
    
    Args:
        state_type (str): Type of regulatory states to return. Options:
            - 'promoter': Promoter-specific states
            - 'enhancer': Enhancer-specific states
            - 'bivalent': Bivalent states
            - 'active': All active states
            - 'all': All regulatory states
    
    Returns:
        dict: Dictionary containing EID lists and regulatory states
    """
    # Define EID categories
    eid_categories = {
        'fetal': [
            "E017", "E093", "E082", "E081", "E089", "E090", "E083", 
            "E092", "E085", "E084", "E086", "E088", "E080"
        ],
        'stem_cell': [
            "E002", "E008", "E001", "E015", "E014", "E016", "E003", "E024",
            "E020", "E019", "E018", "E021", "E022", "E007", "E009", "E010",
            "E013", "E012", "E011", "E004", "E005", "E006"
        ],
        'adult': [
            "E062", "E034", "E045", "E033", "E044", "E043", "E039", "E041",
            "E042", "E040", "E037", "E048", "E038", "E047", "E029", "E031",
            "E035", "E051", "E050", "E036", "E032", "E046", "E030", "E025",
            "E052", "E055", "E056", "E059", "E061", "E057", "E058", "E028",
            "E027", "E054", "E053", "E112", "E071", "E074", "E068", "E069",
            "E072", "E067", "E073", "E070", "E063", "E100", "E108", "E107",
            "E104", "E095", "E105", "E065", "E078", "E076", "E103", "E111",
            "E109", "E106", "E075", "E101", "E102", "E110", "E077", "E079",
            "E094", "E099", "E097", "E087", "E091", "E066", "E098", "E096",
            "E113", "E114", "E119", "E120", "E121", "E122", "E124", "E125",
            "E126", "E127", "E128", "E129"
        ]
    }
    
    # Calculate developmental EIDs
    eid_categories['developmental'] = eid_categories['fetal'] + eid_categories['stem_cell']
    
    # Define different regulatory state configurations
    regulatory_states = {
        'promoter': [
            "1_TssA",      # Active TSS
            "2_TssFlnk",   # Flanking TSS
            "3_TssFlnkU",  # Flanking TSS Upstream
            "4_TssFlnkD",  # Flanking TSS Downstream
        ],
        'enhancer': [
            "7_EnhG1",     # Genic Enhancer 1
            "8_EnhG2",     # Genic Enhancer 2
            "9_EnhA1",     # Active Enhancer 1
            "10_EnhA2",    # Active Enhancer 2
            "11_EnhWk",    # Weak Enhancer
            "15_EnhBiv"    # Bivalent Enhancer
        ],
        'bivalent': [
            "11_EnhWk",    # Weak Enhancer
            "14_TssBiv",   # Bivalent TSS
            "15_EnhBiv"    # Bivalent Enhancer
        ],
        'active': [
            "1_TssA",      # Active TSS
            "2_TssFlnk",   # Flanking TSS
            "3_TssFlnkU",  # Flanking TSS Upstream
            "4_TssFlnkD",  # Flanking TSS Downstream
            "5_Tx",        # Strong Transcription
            "6_TxWk",      # Weak Transcription
            "7_EnhG1",     # Genic Enhancer 1
            "8_EnhG2",     # Genic Enhancer 2
            "9_EnhA1",     # Active Enhancer 1
            "10_EnhA2",    # Active Enhancer 2
            "11_EnhWk",    # Weak Enhancer
            "14_TssBiv",   # Bivalent TSS
            "15_EnhBiv"    # Bivalent Enhancer
        ],
        'all': [
            "1_TssA",      # Active TSS
            "2_TssFlnk",   # Flanking TSS
            "3_TssFlnkU",  # Flanking TSS Upstream
            "4_TssFlnkD",  # Flanking TSS Downstream
            "5_Tx",        # Strong Transcription
            "6_TxWk",      # Weak Transcription
            "7_EnhG1",     # Genic Enhancer 1
            "8_EnhG2",     # Genic Enhancer 2
            "9_EnhA1",     # Active Enhancer 1
            "10_EnhA2",    # Active Enhancer 2
            "11_EnhWk",    # Weak Enhancer
            "12_ReprPC",   # Repressed Polycomb
            "13_ReprPCWk", # Weak Repressed Polycomb
            "14_TssBiv",   # Bivalent TSS
            "15_EnhBiv",   # Bivalent Enhancer
            "16_Quies"     # Quiescent/Low
        ]
    }
    
    # Validate state type
    if state_type not in regulatory_states:
        raise ValueError(f"Invalid state_type: {state_type}. Must be one of: {list(regulatory_states.keys())}")
    
    # Log the configuration
    logger.info(f"Using regulatory state type: {state_type}")
    logger.info(f"Number of states: {len(regulatory_states[state_type])}")
    logger.info(f"States: {', '.join(regulatory_states[state_type])}")
    
    return {
        'eid_categories': eid_categories,
        'regulatory_states': regulatory_states[state_type]
    }


def define_specificity_categories(output_dir, regulatory_regions_dir, mode='cfpeak'):
    """Define specificity categories for regulatory regions based on their activity patterns.
    
    This function analyzes regulatory regions from fetal, stem cell, and adult tissues to identify
    regions with specific activity patterns. It creates BED files for different categories of
    regulatory regions based on their tissue-specific activity.
    
    Categories defined:
    1. Fetal-specific: Regions active only in fetal tissues
    2. Stem cell-specific: Regions active only in stem cells
    3. Fetal-stem shared specific: Regions active in both fetal and stem cells, but not in adult
    4. Developmental-specific: Regions active in developmental stages (fetal and/or stem)
    5. Adult-specific: Regions active only in adult tissues
    6. Fetal-stem shared: Regions active in both fetal and stem cells
    
    Args:
        output_dir (str): Base directory for output files
        regulatory_regions_dir (str): Directory containing the regulatory region BED files
        mode (str): Analysis mode - 'cfpeak' or 'cqtl' (default: 'cfpeak')
    
    Returns:
        dict: Dictionary mapping category names to their BED file paths
    """
    logger.info("Defining specificity categories...")
    
    # Create output directories
    specificity_dir = Path(regulatory_regions_dir) / "specificity_categories"
    specificity_dir.mkdir(exist_ok=True)
    
    # Check if all category files already exist
    expected_categories = [
        "fetal_specific", "stem_cell_specific",
        "developmental_specific", "adult_specific"
    ]
    
    all_files_exist = all(
        (specificity_dir / f"{cat}.bed").exists()
        for cat in expected_categories
    )
    
    if all_files_exist:
        logger.info("All category specificity files already exist. Using existing files.")
        return {cat: specificity_dir / f"{cat}.bed" for cat in expected_categories}
    
    # Get the processed regulatory region files
    logger.info("Step 1/5: Loading regulatory region files...")

    # Look in the nested regulatory_regions directory
    regulatory_bed_dir = Path(regulatory_regions_dir)
    fetal_file = regulatory_bed_dir / "fetal_regulatory_regions.bed"
    stem_file = regulatory_bed_dir / "stem_cell_regulatory_regions.bed"
    adult_file = regulatory_bed_dir / "adult_regulatory_regions.bed"
    
    if not all(f.exists() for f in [fetal_file, stem_file, adult_file]):
        raise FileNotFoundError(f"One or more regulatory region files not found in {regulatory_bed_dir}")
    
    try:
        # Load regulatory regions as BedTools
        logger.info("Step 2/5: Converting files to BedTool format...")
        fetal_bed = BedTool(str(fetal_file))
        stem_bed = BedTool(str(stem_file))
        adult_bed = BedTool(str(adult_file))
        
        # First, create a mapping of coordinates to full records from original file
        logger.info("Step 3/5: Creating coordinate mapping from original file...")
        coord_to_record = {}
        original_file = Path(output_dir) / f"original_{'cqtls' if mode == 'cqtl' else 'cfpeaks'}_with_id.bed"
        
        with open(original_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    coord_key = f"{fields[0]}:{fields[1]}-{fields[2]}"
                    coord_to_record[coord_key] = line.strip()
        
        # Define categories and save with full record information
        logger.info("Step 4/5: Defining and saving specificity categories...")
        categories = {
            "fetal_specific": fetal_bed.intersect(stem_bed, v=True).intersect(adult_bed, v=True),
            "stem_cell_specific": stem_bed.intersect(fetal_bed, v=True).intersect(adult_bed, v=True),
            "developmental_specific": (fetal_bed.cat(stem_bed, postmerge=False).sort().merge()).intersect(adult_bed, v=True),
            "adult_specific": adult_bed.intersect(fetal_bed, v=True).intersect(stem_bed, v=True)
        }
        
        # Calculate total regions from original files
        total_regions = len(fetal_bed) + len(stem_bed) + len(adult_bed)
        logger.info(f"Total regions from original files: {total_regions}")
        
        # Save each category with full record information
        category_files = {}
        for category_name, bed in categories.items():
            # Save basic BED file for internal use
            output_file = specificity_dir / f"{category_name}.bed"
            bed.saveas(str(output_file))
            category_files[category_name] = output_file
            
            count = len(bed)
            percentage = (count / total_regions * 100) if total_regions > 0 else 0
            logger.info(f"\n{category_name.replace('_', ' ').title()}:")
            logger.info(f"  Count: {count}")
            logger.info(f"  Percentage: {percentage:.1f}%")
            logger.info(f"  Basic BED saved to: {output_file}")
        
        return category_files
        
    except Exception as e:
        logger.error(f"Error defining specificity categories: {e}")
        raise

def print_category_eids(bed_file, category_name):
    """Extract and print EIDs from a BED file for a specific category.
    
    Args:
        bed_file (str): Path to the BED file
        category_name (str): Name of the category for logging
    """
    try:
        eids = set()
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('track'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    eids.add(fields[3])
        
        eids = sorted(list(eids))
        logger.info(f"\n{category_name} EIDs ({len(eids)} total):")
        logger.info("-" * 50)
        logger.info(", ".join(eids))
        
        # Save EIDs to a separate file
        eid_file = Path(bed_file).parent / f"{Path(bed_file).stem}_eids.txt"
        with open(eid_file, 'w') as f:
            f.write("\n".join(eids))
        logger.info(f"EIDs saved to: {eid_file}")
        
    except Exception as e:
        logger.error(f"Error extracting EIDs from {bed_file}: {e}")


def process_cfpeak_file(cfpeak_file, output_dir):
    """Process cfPeak file and convert to BED format with identifiers."""
    logger.info(f"Processing cfPeak file: {cfpeak_file}")
    
    try:
        # Create output directory
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Read and check format
        with open(cfpeak_file, 'r') as f:
            first_line = f.readline().strip()
        
        if not first_line:
            raise ValueError("Empty cfPeak file")
        
        # Create a copy of the original file with identifiers
        original_with_id = output_dir / "original_cfpeaks_with_id.bed"
        
        with open(cfpeak_file, 'r') as infile, \
             open(original_with_id, 'w') as outfile:
            for i, line in enumerate(infile, 1):
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    try:
                        chrom = fields[0]
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        start = int(fields[1])
                        end = int(fields[2])
                        if start >= end:
                            logger.warning(f"Invalid coordinates in line: {line.strip()}")
                            continue
                        
                        # Create identifier: cfpeak_<index>
                        identifier = f"cfpeak_{i}"
                        
                        # Write line with identifier
                        output_fields = [chrom, str(start), str(end), identifier]
                        
                        # Add any additional fields from the original file
                        if len(fields) > 3:
                            output_fields.extend(fields[3:])
                        
                        outfile.write('\t'.join(output_fields) + '\n')
                        
                    except ValueError as e:
                        logger.warning(f"Invalid line format: {line.strip()} - {e}")
                        continue
        
        # Save the original file without identifiers for reference
        original_copy = output_dir / "original_cfpeaks.bed"
        with open(cfpeak_file, 'r') as infile, \
             open(original_copy, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    try:
                        chrom = fields[0]
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        start = int(fields[1])
                        end = int(fields[2])
                        if start >= end:
                            continue
                        output_fields = [chrom, str(start), str(end)]
                        if len(fields) > 3:
                            output_fields.extend(fields[3:])
                        outfile.write('\t'.join(output_fields) + '\n')
                    except ValueError:
                        continue
        
        # Convert to BedTool using the file with identifiers
        cfpeaks_bed = BedTool(str(original_with_id))
        total_cfpeaks = len(cfpeaks_bed)
        
        if total_cfpeaks == 0:
            raise ValueError("No valid peaks found in cfPeak file")
        
        logger.info(f"Total cfPeaks: {total_cfpeaks}")
        logger.info(f"Original file with identifiers saved to: {original_with_id}")
        logger.info(f"Original file copy saved to: {original_copy}")
        
        return cfpeaks_bed, total_cfpeaks
        
    except Exception as e:
        logger.error(f"Error processing cfPeak file: {e}")
        raise

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

def process_custom_random_bg(custom_bg_file, base_dir=None, target_count=None, seed=50, pval_threshold=None, qval_threshold=None, mode='bed', perform_permutation=False, n_permutations=100, return_full_bedtool=False):
    """Process custom random background file and convert to BED format.
    
    Args:
        custom_bg_file (str): Path to custom random background file
        base_dir (str): Base directory for output files
        target_count (int): Target number of regions to sample
        seed (int): Random seed for reproducibility (Default: 50)
        pval_threshold (float): P-value threshold for filtering (default: 0.5)
        qval_threshold (float): Q-value threshold for filtering (default: None)
        mode (str): File type mode. Options:
            - 'bed': Standard BED format file (default)
            - 'cqtl': cQTL format file
            - 'cfpeak': cfPeak format file
        perform_permutation (bool): Whether to perform permutation sampling
        n_permutations (int): Number of permutations to perform
        return_full_bedtool (bool): If True, also return the full (unsampled) BedTool object.
    
    Returns:
        tuple: (BedTool object, dict) containing the processed random background and metadata including seed information
    """
    logger.info(f"Processing custom random background file: {custom_bg_file}")
    
    # Create random background directory
    random_bg_dir = Path(base_dir) / "random_background"
    random_bg_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize metadata dictionary
    metadata = {
        'seed': seed,
        'target_count': target_count,
        'mode': mode,
        'pval_threshold': pval_threshold,
        'qval_threshold': qval_threshold,
        'perform_permutation': perform_permutation,
        'n_permutations': n_permutations,
        'filtering_method': 'q-value' if qval_threshold is not None else 'p-value'
    }
    
    try:
        # First try reading with header to check column names
        try:
            logger.info(f"Attempting to read file with headers: {custom_bg_file}")
            df = pd.read_csv(custom_bg_file, sep='\t', low_memory=False)
            has_header = True
            logger.info(f"Successfully read with headers. Columns: {df.columns.tolist()}")
        except Exception as e:
            logger.info(f"Failed to read with headers: {e}. Trying without headers...")
            # If that fails, read without header
            df = pd.read_csv(custom_bg_file, sep='\t', header=None, low_memory=False)
            has_header = False
            logger.info(f"Successfully read without headers. Number of columns: {len(df.columns)}")
        
        # Process based on mode
        if mode == 'cqtl':
            logger.info("Processing as cQTL file")
            
            if has_header:
                # Check if we have the required columns
                if 'CHR' in df.columns and 'POS' in df.columns:
                    logger.info("Found CHR and POS columns in header")
                    # If p-value or q-value threshold is provided, filter by significance
                    original_count = len(df)
                    
                    # Log the thresholds being used
                    logger.info("\nThreshold Information:")
                    logger.info(f"  P-value threshold: {pval_threshold}")
                    logger.info(f"  Q-value threshold: {qval_threshold}")
                    
                    if qval_threshold is not None and 'comb.q' in df.columns:
                        logger.info("\nUsing Q-value based filtering (takes precedence over p-value)")
                        logger.info(f"Filtering cQTLs by q-value threshold {qval_threshold}")
                        df = df[df['comb.q'] <= qval_threshold]
                        filtered_count = len(df)
                        logger.info(f"Filtered cQTLs by q-value threshold {qval_threshold}:")
                        logger.info(f"  Original count: {original_count}")
                        logger.info(f"  After filtering: {filtered_count}")
                        logger.info(f"  Removed: {original_count - filtered_count} cQTLs")
                    elif pval_threshold is not None and 'pval' in df.columns:
                        logger.info("\nUsing P-value based filtering")
                        logger.info(f"Filtering cQTLs by p-value threshold {pval_threshold}")
                        df = df[df['pval'] <= pval_threshold]
                        filtered_count = len(df)
                        logger.info(f"Filtered cQTLs by p-value threshold {pval_threshold}:")
                        logger.info(f"  Original count: {original_count}")
                        logger.info(f"  After filtering: {filtered_count}")
                        logger.info(f"  Removed: {original_count - filtered_count} cQTLs")
                    
                    # Create BED format DataFrame from filtered data
                    bed_df = pd.DataFrame()
                    bed_df['chrom'] = 'chr' + df['CHR'].astype(str)
                    bed_df['start'] = df['POS'] - 1  # BED is 0-based
                    bed_df['end'] = df['POS']        # End is exclusive in BED format
                else:
                    # If no CHR/POS columns, assume first two columns are chromosome and position
                    logger.info("No CHR/POS columns found, using first two columns")
                    bed_df = pd.DataFrame()
                    bed_df['chrom'] = 'chr' + df.iloc[:, 0].astype(str)
                    bed_df['start'] = df.iloc[:, 1] - 1
                    bed_df['end'] = df.iloc[:, 1]
            else:
                # No headers, assume first two columns are chromosome and position
                logger.info("No headers, using first two columns")
                bed_df = pd.DataFrame()
                bed_df['chrom'] = 'chr' + df.iloc[:, 0].astype(str)
                bed_df['start'] = df.iloc[:, 1] - 1
                bed_df['end'] = df.iloc[:, 1]
            
            # Add identifier column
            bed_df['identifier'] = [f"random_{i+1}" for i in range(len(bed_df))]
            
            # Save processed BED file with identifiers
            output_file = random_bg_dir / "processed_random.bed"
            bed_df.to_csv(output_file, sep='\t', header=False, index=False)
            
            logger.info(f"Converted cQTL to BED format and saved to: {output_file}")
            logger.info(f"Final bed_df shape: {bed_df.shape}")
            
            # Load as BedTool (full set)
            full_bed = BedTool(str(output_file))
            
        elif mode == 'cfpeak':
            logger.info("Processing as cfPeak file")
            # Create BED format DataFrame
            bed_df = pd.DataFrame()
            
            if has_header:
                # Check if we have the required columns
                if 'CHR' in df.columns and 'START' in df.columns and 'END' in df.columns:
                    logger.info("Found CHR, START, and END columns in header")
                    bed_df['chrom'] = 'chr' + df['CHR'].astype(str)
                    bed_df['start'] = df['START']
                    bed_df['end'] = df['END']
                else:
                    # If no standard columns, assume first three columns are chromosome, start, end
                    logger.info("No standard columns found, using first three columns")
                    bed_df['chrom'] = 'chr' + df.iloc[:, 0].astype(str)
                    bed_df['start'] = df.iloc[:, 1]
                    bed_df['end'] = df.iloc[:, 2]
            else:
                # No headers, assume first three columns are chromosome, start, end
                logger.info("No headers, using first three columns")
                bed_df['chrom'] = 'chr' + df.iloc[:, 0].astype(str)
                bed_df['start'] = df.iloc[:, 1]
                bed_df['end'] = df.iloc[:, 2]
            
            # Add identifier column
            bed_df['identifier'] = [f"random_{i+1}" for i in range(len(bed_df))]
            
            # Save processed BED file with identifiers
            output_file = random_bg_dir / "processed_random.bed"
            bed_df.to_csv(output_file, sep='\t', header=False, index=False)
            
            logger.info(f"Converted cfPeak to BED format and saved to: {output_file}")
            logger.info(f"Filtered bed_df shape: {bed_df.shape}")
            
        else:  # mode == 'bed'
            logger.info("Processing as BED file")
            # Ensure chromosome format
            if has_header:
                chrom_col = df.columns[0]
                df[chrom_col] = df[chrom_col].astype(str)
                df[chrom_col] = df[chrom_col].apply(lambda x: x if x.startswith('chr') else 'chr' + x)
            else:
                df[0] = df[0].astype(str)
                df[0] = df[0].apply(lambda x: x if x.startswith('chr') else 'chr' + x)
            
            # Add identifier column
            df['identifier'] = [f"random_{i+1}" for i in range(len(df))]
            
            # Save processed BED file with identifiers
            output_file = random_bg_dir / "processed_random.bed"
            if has_header:
                df.iloc[:, :3].join(df['identifier']).to_csv(output_file, sep='\t', header=False, index=False)
            else:
                df.iloc[:, :3].join(df['identifier']).to_csv(output_file, sep='\t', header=False, index=False)
            
            logger.info(f"Processed BED file saved to: {output_file}")
            logger.info(f"Final df shape: {df.shape}")
        
        # Sample to target count if specified
        if target_count is not None:
            # Load as BedTool
            bed_tool = BedTool(str(output_file))
            original_count = len(bed_tool)
            
            if target_count < original_count:
                logger.info(f"Randomly sampling {target_count} items from {original_count} total items...")
                
                # Get all intervals
                intervals = list(bed_tool)
                
                # Set random seed for reproducibility
                random.seed(seed)
                
                if perform_permutation:
                    # Sample multiple times for permutation test
                    all_sampled_intervals = []
                    for _ in range(n_permutations):
                        sampled_indices = random.sample(range(len(intervals)), target_count)
                        all_sampled_intervals.append([intervals[i] for i in sampled_indices])
                    # Use first permutation as the main sample
                    sampled_intervals = all_sampled_intervals[0]
                    
                    # Save all permutations to separate files
                    for i, perm_intervals in enumerate(all_sampled_intervals):
                        perm_file = random_bg_dir / f"permutation_{i+1}.bed"
                        with open(perm_file, 'w') as f:
                            for interval in perm_intervals:
                                f.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")
                        logger.info(f"Saved permutation {i+1} to: {perm_file}")
                else:
                    # Single random sampling
                    sampled_indices = random.sample(range(len(intervals)), target_count)
                    sampled_intervals = [intervals[i] for i in sampled_indices]
                
                # Save to temporary file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
                    for interval in sampled_intervals:
                        tmp.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")
                    tmp_path = tmp.name
                
                # Update BedTool with sampled intervals
                bed_tool = BedTool(tmp_path)
                
                # Save sampled intervals to a new file for reference
                sampled_file = random_bg_dir / "sampled_random.bed"
                bed_tool.saveas(str(sampled_file))
                
                # Create sampled file with identifiers
                sampled_with_id = random_bg_dir / "sampled_random_with_id.bed"
                with open(sampled_with_id, 'w') as f:
                    for i, interval in enumerate(sampled_intervals, 1):
                        f.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\trandom_{i}\n")
                
                logger.info(f"Sampled intervals saved to: {sampled_file}")
                logger.info(f"Sampled intervals with identifiers saved to: {sampled_with_id}")
                
                output_file = sampled_file
        
        # Load final BedTool
        final_bed = BedTool(str(output_file))
        
        # Update metadata with final counts
        metadata.update({
            'original_count': original_count if target_count is not None else len(final_bed),
            'final_count': len(final_bed),
            'output_file': str(output_file)
        })
        
        logger.info(f"Loaded custom random background:")
        logger.info(f"  Original count: {metadata['original_count']}")
        logger.info(f"  Final count: {metadata['final_count']}")
        logger.info(f"  Target count: {target_count if target_count else 'All items'}")
        logger.info(f"  Random seed used: {seed}")
        if perform_permutation:
            logger.info(f"  Performed {n_permutations} permutations")
        
        if return_full_bedtool:
            return final_bed, metadata, full_bed
        else:
            return final_bed, metadata
        
    except Exception as e:
        logger.error(f"Error processing custom random background file: {e}")
        raise

def perform_overlap_analysis(cfpeak_cqtl_bed, custom_regions_bed, custom_random_bg, output_dir, category_name=None, mode='bed', source_files=None, perform_permutation=False, n_permutations=100):
    """Perform basic overlap analysis between input regions and custom regions."""
    logger.info(f"Performing overlap analysis for {mode} mode...")
    
    try:
        # Unpack custom_random_bg tuple
        if isinstance(custom_random_bg, tuple) and len(custom_random_bg) == 3:
            random_bg_bed, random_bg_metadata, full_random_bed = custom_random_bg
        else:
            random_bg_bed, random_bg_metadata = custom_random_bg
            full_random_bed = random_bg_bed
        
        # Calculate overlaps using bedtools unique
        logger.info(f"Calculating overlaps for input regions with {category_name} regions...")
        overlaps = cfpeak_cqtl_bed.intersect(custom_regions_bed, u=True)
        random_overlaps = random_bg_bed.intersect(custom_regions_bed, u=True)
        
        # Get counts
        n_overlaps = len(overlaps)
        n_random_overlaps = len(random_overlaps)
        total_regions = len(cfpeak_cqtl_bed)
        total_random = len(random_bg_bed)
        total_custom_regions = len(custom_regions_bed)
        
        logger.info(f"Input regions: {total_regions} total, {n_overlaps} overlaps")
        logger.info(f"Random background: {total_random} total, {n_random_overlaps} overlaps")
        
        # Calculate overlap fraction
        overlap_fraction = n_overlaps / total_regions if total_regions > 0 else 0
        
        # Perform Fisher's exact test with confidence intervals
        contingency_table = np.array([
            [n_overlaps, total_regions - n_overlaps],
            [n_random_overlaps, total_random - n_random_overlaps]
        ])
        
        # Calculate odds ratio and p-value
        odds_ratio, p_value = stats.fisher_exact(contingency_table)
        
        # Calculate 95% confidence intervals for odds ratio
        alpha = 0.05
        z = stats.norm.ppf(1 - alpha/2)
        log_or = np.log(odds_ratio)
        se_log_or = np.sqrt(sum(1/x for x in contingency_table.flatten()))
        ci_lower = np.exp(log_or - z * se_log_or)
        ci_upper = np.exp(log_or + z * se_log_or)
        
        # Calculate fold enrichment
        fold_enrichment = (n_overlaps / total_regions) / (n_random_overlaps / total_random) if n_random_overlaps > 0 else 0
        
        # Create results dictionary for detailed summary
        results_dict = {
            'n_overlaps': n_overlaps,
            'total_regions': total_regions,
            'n_random_overlaps': n_random_overlaps,  # Add this key
            'total_random': total_random,
            'total_custom_regions': total_custom_regions,
            'overlap_fraction': overlap_fraction,
            'fold_enrichment': fold_enrichment,
            'fisher_p_value': p_value,
            'odds_ratio': odds_ratio,
            'fishers_exact': {
                'p_value': p_value,
                'odds_ratio': odds_ratio,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper
            }
        }
        
        # Perform permutation analysis if requested
        permutation_results = None
        if perform_permutation and random_bg_metadata.get('perform_permutation', False):
            logger.info(f"Performing {n_permutations} permutations...")
            permutation_fold_enrichments = []
            permutation_contingency_tables = []
            permutation_odds_ratios = []
            permutation_odds_ratio_cis = []
            permutation_p_values = []  # Add list to store p-values
            
            # Get the random background directory
            random_bg_dir = Path(output_dir) / "random_background"
            
            # Process each permutation file
            for i in range(1, n_permutations + 1):
                perm_file = random_bg_dir / f"permutation_{i}.bed"
                if not perm_file.exists():
                    logger.warning(f"Permutation file {perm_file} not found. Skipping...")
                    continue
                
                logger.info(f"\nProcessing permutation {i}:")
                logger.info(f"Loading permutation file: {perm_file}")
                
                # Load permutation as BedTool
                perm_bed = BedTool(str(perm_file))
                total_perm = len(perm_bed)
                
                # Calculate overlaps for this permutation
                perm_overlaps = perm_bed.intersect(custom_regions_bed, u=True)
                n_perm_overlaps = len(perm_overlaps)
                
                logger.info(f"Permutation {i} regions: {total_perm} total, {n_perm_overlaps} overlaps")
                
                # Calculate fold enrichment for this permutation
                perm_fold_enrichment = (n_overlaps / total_regions) / (n_perm_overlaps / total_perm) if n_perm_overlaps > 0 else 0
                permutation_fold_enrichments.append(perm_fold_enrichment)
                
                # Store contingency table for this permutation
                perm_contingency = np.array([
                    [n_overlaps, total_regions - n_overlaps],
                    [n_perm_overlaps, total_perm - n_perm_overlaps]
                ])
                permutation_contingency_tables.append(perm_contingency)
                
                # Calculate Fisher's exact test for this permutation
                perm_odds_ratio, perm_p_value = stats.fisher_exact(perm_contingency)
                permutation_p_values.append(perm_p_value)  # Store p-value
                
                # Calculate confidence intervals for odds ratio
                alpha = 0.05
                z = stats.norm.ppf(1 - alpha/2)
                log_or = np.log(perm_odds_ratio)
                se_log_or = np.sqrt(sum(1/x for x in perm_contingency.flatten()))
                ci_lower = np.exp(log_or - z * se_log_or)
                ci_upper = np.exp(log_or + z * se_log_or)
                
                permutation_odds_ratios.append(perm_odds_ratio)
                permutation_odds_ratio_cis.append((ci_lower, ci_upper))
                
                logger.info(f"Permutation {i} fold enrichment: {perm_fold_enrichment:.4f}")
                logger.info(f"Permutation {i} odds ratio: {perm_odds_ratio:.4f} (95% CI: {ci_lower:.4f}-{ci_upper:.4f})")
                logger.info(f"Permutation {i} p-value: {perm_p_value:.4e}")
                logger.info(f"Total original overlaps: {n_overlaps}")
            
            # Calculate statistics for permutation results
            mean_permutation_fold = np.mean(permutation_fold_enrichments)
            std_permutation_fold = np.std(permutation_fold_enrichments)
            
            # Filter out infinite values from odds ratios
            finite_permutation_odds = [x for x in permutation_odds_ratios if np.isfinite(x)]
            if not finite_permutation_odds:
                logger.warning("No finite odds ratios found in permutations. Skipping odds ratio statistics.")
                mean_permutation_odds = np.nan
                std_permutation_odds = np.nan
            else:
                mean_permutation_odds = np.mean(finite_permutation_odds)
                std_permutation_odds = np.std(finite_permutation_odds)
            
            permutation_p_value = (sum(1 for x in permutation_fold_enrichments if x >= fold_enrichment) + 1) / (n_permutations + 1)
            
            # Create plots
            plots_dir = Path(output_dir) / "plots"
            plots_dir.mkdir(exist_ok=True)
            
            # Set style for publication-quality plots
            try:
                plt.style.use('seaborn-v0_8-whitegrid')
            except OSError:
                # Fallback to a simpler style if seaborn is not available
                plt.style.use('default')
                # Manually set some style parameters for clean look
                plt.rcParams.update({
                    'axes.grid': True,
                    'grid.alpha': 0.3,
                    'axes.facecolor': 'white',
                    'figure.facecolor': 'white',
                    'axes.edgecolor': 'black',
                    'axes.linewidth': 1.0
                })
            
            try:
                # 2. Histogram plots
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
                
                # Fold enrichment histogram
                ax1.hist(permutation_fold_enrichments, bins=20, alpha=0.7, color='#3498db')
                ax1.axvline(mean_permutation_fold, color='red', linestyle='--', 
                           label=f'Mean: {mean_permutation_fold:.2f}')
                ax1.axvline(fold_enrichment, color='green', linestyle='--', 
                           label=f'Observed: {fold_enrichment:.2f}')
                ax1.set_title('Permutation Fold Enrichments')
                ax1.set_xlabel('Fold Enrichment')
                ax1.set_ylabel('Count')
                ax1.legend()
                
                # Odds ratio histogram (only if we have finite values)
                if finite_permutation_odds:
                    # Calculate histogram bins based on finite values
                    min_val = min(finite_permutation_odds)
                    max_val = max(finite_permutation_odds)
                    bins = np.linspace(min_val, max_val, 21)  # 20 bins + 1 edge
                    
                    ax2.hist(finite_permutation_odds, bins=bins, alpha=0.7, color='#3498db')
                    ax2.axvline(mean_permutation_odds, color='red', linestyle='--', 
                               label=f'Mean: {mean_permutation_odds:.2f}')
                    ax2.axvline(odds_ratio, color='green', linestyle='--', 
                               label=f'Observed: {odds_ratio:.2f}')
                    ax2.set_title('Permutation Odds Ratios')
                    ax2.set_xlabel('Odds Ratio')
                    ax2.set_ylabel('Count')
                    ax2.legend()
                else:
                    ax2.text(0.5, 0.5, 'No finite odds ratios available', 
                            ha='center', va='center', transform=ax2.transAxes)
                    ax2.set_title('Permutation Odds Ratios')
                
                plt.tight_layout()
                plt.savefig(plots_dir / f"{category_name}_histogram_plot.pdf", dpi=300, bbox_inches='tight')
                plt.close()
            except Exception as e:
                logger.error(f"Error creating histograms: {e}")
                logger.warning("Skipping histogram creation due to error")
            
            permutation_results = {
                'mean_fold_enrichment': mean_permutation_fold,
                'std_fold_enrichment': std_permutation_fold,
                'mean_odds_ratio': mean_permutation_odds,
                'std_odds_ratio': std_permutation_odds,
                'p_value': permutation_p_value,
                'all_fold_enrichments': permutation_fold_enrichments,
                'all_odds_ratios': permutation_odds_ratios,
                'all_odds_ratio_cis': permutation_odds_ratio_cis,
                'all_p_values': permutation_p_values,
                'contingency_tables': permutation_contingency_tables,
                'base_seed': random_bg_metadata['seed'],
                'n_permutations': n_permutations,
                'n_finite_odds_ratios': len(finite_permutation_odds)
            }
        
        # Create detailed statistics dictionary
        stats_dict = OrderedDict([
            ('category', category_name if category_name else 'unknown'),
            ('total_peaks', total_regions),
            ('overlapping_peaks', n_overlaps), 
            ('percentage', (n_overlaps / total_regions * 100) if total_regions > 0 else 0),
            ('fold_enrichment', fold_enrichment),
            ('total_regions', total_custom_regions),
            ('fisher_p_value', p_value),
            ('odds_ratio', odds_ratio),
            ('odds_ratio_ci_lower', ci_lower),
            ('odds_ratio_ci_upper', ci_upper)
        ])
        
        # Add permutation results if available
        if permutation_results:
            stats_dict.update({
                'permutation_mean_fold': permutation_results['mean_fold_enrichment'],
                'permutation_std_fold': permutation_results['std_fold_enrichment'],
                'permutation_mean_odds': permutation_results['mean_odds_ratio'],
                'permutation_std_odds': permutation_results['std_odds_ratio'],
                'permutation_p_value': permutation_results['p_value'],
                'permutation_base_seed': permutation_results['base_seed']
            })
        
        # Save detailed statistics to CSV
        output_dir = Path(output_dir)
        stats_file = output_dir / 'detailed_statistics.csv'
        
        # Convert stats to DataFrame and save
        stats_df = pd.DataFrame([stats_dict])
        
        # If file exists, append to it, otherwise create new
        if stats_file.exists():
            existing_df = pd.read_csv(stats_file)
            stats_df = pd.concat([existing_df, stats_df], ignore_index=True)
            
            # Ensure category is the first column by reordering
            column_order = ['category'] + [col for col in stats_df.columns if col != 'category']
            stats_df = stats_df[column_order]
        stats_df.to_csv(stats_file, index=False)
        
        # Add permutation results if available
        if permutation_results:
            results_dict['permutation_analysis'] = permutation_results
        
        # Create detailed summary file
        create_detailed_summary(results_dict, output_dir, category_name, source_files)
        
        # Print results
        print(f"\nResults for category: {category_name}")
        print("=" * 50)
        print(f"Total peaks: {total_regions}")
        print(f"Overlapping peaks: {n_overlaps}")
        print(f"Percentage: {stats_dict['percentage']:.2f}%")
        print(f"Fold enrichment: {fold_enrichment:.2f}")
        print(f"Fisher's exact test p-value: {p_value:.2e}")
        print(f"Odds ratio: {odds_ratio:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f})")
        
        if permutation_results:
            print("\nPermutation Analysis Results:")
            print(f"Mean fold enrichment: {permutation_results['mean_fold_enrichment']:.2f} ± {permutation_results['std_fold_enrichment']:.2f}")
            print(f"Mean odds ratio: {permutation_results['mean_odds_ratio']:.4f} ± {permutation_results['std_odds_ratio']:.4f}")
            print(f"Permutation p-value: {permutation_results['p_value']:.4e}")
            print(f"Base seed used: {permutation_results['base_seed']}")
            print(f"Number of permutations: {permutation_results['n_permutations']}")
        
        print("=" * 50)
        
        logger.info(f"Detailed statistics saved to: {stats_file}")
        
        # --- NEW: Save full random background Fisher's test ---
        # Use full_random_bed for the full random background stats
        save_full_random_fisher_stats(
            cfpeak_cqtl_bed,
            full_random_bed,
            custom_regions_bed,
            output_dir,
            category_name=category_name
        )
        
        return results_dict
        
    except Exception as e:
        logger.error(f"Error performing overlap analysis: {e}")
        raise

def save_full_random_fisher_stats(input_bed, full_random_bed, custom_regions_bed, output_dir, category_name=None):
    """Save Fisher's test results using the full set of random cQTLs as reference."""
    from scipy import stats
    import numpy as np
    import pandas as pd
    from collections import OrderedDict
    from pathlib import Path

    # Calculate overlaps
    n_input = len(input_bed)
    n_full_random = len(full_random_bed)
    n_input_overlap = len(input_bed.intersect(custom_regions_bed, u=True))
    n_full_random_overlap = len(full_random_bed.intersect(custom_regions_bed, u=True))
    n_total_regions = len(custom_regions_bed)

    # Fisher's exact test
    table = np.array([
        [n_input_overlap, n_input - n_input_overlap],
        [n_full_random_overlap, n_full_random - n_full_random_overlap]
    ])
    odds_ratio, p_value = stats.fisher_exact(table)
    # Confidence interval
    alpha = 0.05
    z = stats.norm.ppf(1 - alpha/2)
    log_or = np.log(odds_ratio) if odds_ratio > 0 else 0
    se_log_or = np.sqrt(np.sum(1 / table))
    ci_lower = np.exp(log_or - z * se_log_or)
    ci_upper = np.exp(log_or + z * se_log_or)
    # Fold enrichment
    fold_enrichment = (n_input_overlap / n_input) / (n_full_random_overlap / n_full_random) if n_full_random_overlap > 0 and n_full_random > 0 else float('inf')

    # Calculate percentage of overlapping input peaks
    percentage = (n_input_overlap / n_input * 100) if n_input > 0 else 0

    # Prepare output
    stats_dict = OrderedDict([
        ('category', category_name if category_name else 'unknown'),
        ('total_peaks', n_input),
        ('overlapping_peaks', n_input_overlap),
        ('percentage', percentage),
        ('full_random_total', n_full_random),
        ('full_random_overlap', n_full_random_overlap),
        ('fold_enrichment', fold_enrichment),
        ('total_regions', n_total_regions),
        ('odds_ratio', odds_ratio),
        ('odds_ratio_ci_lower', ci_lower),
        ('odds_ratio_ci_upper', ci_upper),
        ('fisher_p_value', p_value)
    ])
    output_dir = Path(output_dir)
    stats_file = output_dir / 'detailed_statistics_full_random.csv'
    stats_df = pd.DataFrame([stats_dict])
    # If file exists, append to it, otherwise create new
    if stats_file.exists():
        existing_df = pd.read_csv(stats_file)
        stats_df = pd.concat([existing_df, stats_df], ignore_index=True)
        # Ensure category is the first column by reordering
        column_order = ['category'] + [col for col in stats_df.columns if col != 'category']
        stats_df = stats_df[column_order]
    stats_df.to_csv(stats_file, index=False)

def create_detailed_summary(results_dict, output_dir, category_name, source_files=None):
    """Create a detailed summary file with contingency tables and source information."""
    output_dir = Path(output_dir)
    summary_file = output_dir / "detailed_analysis_summary.txt"
    
    with open(summary_file, 'a') as f:
        f.write(f"\n{'='*80}\n")
        f.write(f"Analysis Summary for Category: {category_name}\n")
        f.write(f"{'='*80}\n\n")
        
        # Write source file information if available
        if source_files:
            f.write("Source Files:\n")
            f.write("-" * 50 + "\n")
            f.write(f"Input File: {source_files['input_file']}\n")
            f.write(f"Random Background File: {source_files['random_bg_file']}\n")
            f.write(f"Custom Regions File: {source_files['custom_regions_file']}\n\n")
        
        # Write threshold information
        f.write("Threshold Information:\n")
        f.write("-" * 50 + "\n")
        if 'metadata' in results_dict:
            metadata = results_dict['metadata']
            f.write(f"P-value threshold: {metadata.get('pval_threshold', 'Not used')}\n")
            f.write(f"Q-value threshold: {metadata.get('qval_threshold', 'Not used')}\n")
            f.write(f"Filtering method used: {metadata.get('filtering_method', 'Not specified')}\n")
        f.write("\n")
        
        # Write main contingency table
        f.write("Main Contingency Table:\n")
        f.write("-" * 50 + "\n")
        f.write("                    | Overlapping | Non-overlapping | Total\n")
        f.write("  -" * 25 + "\n")
        f.write(f"  Input Regions      | {results_dict['n_overlaps']:11d} | {results_dict['total_regions'] - results_dict['n_overlaps']:14d} | {results_dict['total_regions']:5d}\n")
        f.write(f"  Random Background  | {results_dict['n_random_overlaps']:11d} | {results_dict['total_random'] - results_dict['n_random_overlaps']:14d} | {results_dict['total_random']:5d}\n")
        f.write("  -" * 25 + "\n\n")
        
        # Write statistical results
        f.write("Statistical Results:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Fisher's Exact Test:\n")
        f.write(f"  P-value: {results_dict['fishers_exact']['p_value']:.4e}\n")
        f.write(f"  Odds Ratio: {results_dict['fishers_exact']['odds_ratio']:.4f}\n")
        f.write(f"  95% CI: [{results_dict['fishers_exact']['ci_lower']:.4f}, {results_dict['fishers_exact']['ci_upper']:.4f}]\n\n")
        
        f.write(f"Fold Enrichment: {results_dict['fold_enrichment']:.4f}\n")
        f.write(f"Fisher's Exact Test P-value: {results_dict['fisher_p_value']:.4e}\n")
        
        # Write permutation results if available
        if 'permutation_analysis' in results_dict:
            f.write("Permutation Analysis Results:\n")
            f.write("-" * 50 + "\n")
            f.write(f"Number of Permutations: {len(results_dict['permutation_analysis']['all_fold_enrichments'])}\n")
            f.write(f"Mean Fold Enrichment: {results_dict['permutation_analysis']['mean_fold_enrichment']:.4f} ± {results_dict['permutation_analysis']['std_fold_enrichment']:.4f}\n")
            f.write(f"Mean Odds Ratio: {results_dict['permutation_analysis']['mean_odds_ratio']:.4f} ± {results_dict['permutation_analysis']['std_odds_ratio']:.4f}\n")
            f.write(f"Permutation P-value: {results_dict['permutation_analysis']['p_value']:.4e}\n\n")
            
            # Write individual permutation results
            f.write("Individual Permutation Results:\n")
            for i, (fold, odds, odds_ci, p_val, contingency) in enumerate(zip(
                results_dict['permutation_analysis']['all_fold_enrichments'],
                results_dict['permutation_analysis']['all_odds_ratios'],
                results_dict['permutation_analysis']['all_odds_ratio_cis'],
                results_dict['permutation_analysis']['all_p_values'],
                results_dict['permutation_analysis']['contingency_tables']
            ), 1):
                f.write(f"\nPermutation {i}:\n")
                f.write(f"  Fold Enrichment: {fold:.4f}\n")
                f.write(f"  Odds Ratio: {odds:.4f} (95% CI: {odds_ci[0]:.4f}-{odds_ci[1]:.4f})\n")
                f.write(f"  Fisher's Exact Test P-value: {p_val:.4e}\n")
                f.write("  Contingency Table:\n")
                f.write("                    | Overlapping | Non-overlapping | Total\n")
                f.write("  -" * 25 + "\n")
                f.write(f"  Input Regions      | {contingency[0][0]:11d} | {contingency[0][1]:14d} | {sum(contingency[0]):5d}\n")
                f.write(f"  Permutation {i:3d}    | {contingency[1][0]:11d} | {contingency[1][1]:14d} | {sum(contingency[1]):5d}\n")
                f.write("  -" * 25 + "\n")
            f.write("\n")
        
        # Write source information
        f.write("Source Information:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total Input Regions: {results_dict['total_regions']}\n")
        f.write(f"Total Random Background Regions: {results_dict['total_random']}\n")
        f.write(f"Total Custom Regions: {results_dict['total_custom_regions']}\n\n")
        
        f.write(f"{'='*80}\n\n")

def process_chromhmm_files(chromhmm_dir, category_name, output_dir, reg_states):
    """Process chromHMM files for a specific category.
    
    Args:
        chromhmm_dir (str): Directory containing chromHMM files
        category_name (str): Name of the category (e.g., 'fetal', 'stem_cell', 'adult')
        output_dir (str): Directory to save processed files
        reg_states (dict): Dictionary containing regulatory states and EID categories
    
    Returns:
        str: Path to the processed BED file
    """
    logger.info(f"Processing chromHMM files for category: {category_name}")

    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate regulatory regions dir
        regulatory_bed_dir = Path(output_dir) / "regulatory_regions"
        regulatory_bed_dir.mkdir(exist_ok=True)  # Create directory if it doesn't exist
        output_file = regulatory_bed_dir / f"{category_name}_regulatory_regions.bed"
        
        # Get regulatory states and EID categories from the passed parameter
        states = reg_states['regulatory_states']
        eid_categories = reg_states['eid_categories']
        
        # category_name = "developmental"
        # Get EIDs for this category
        category_eids = eid_categories.get(category_name, [])
        
        if not category_eids:
            logger.warning(f"No EIDs found for {category_name} category")
            return None
        
        print("")
        print("")
        logger.info(f"Regulatory states to filter: {states}")
        logger.info(f"Found {len(category_eids)} EIDs for {category_name} category")
        print("")
        print("")

        all_regions = []

        # For each EID, look for a chromHMM BED file and process it
        for eid in category_eids:
            # Match file pattern: EID_*_core_K27ac_dense.bed or .bed.gz
            bed_file = None
            for ext in [".bed.gz", ".bed"]:
                pattern = f"{eid}_18_core_K27ac_dense{ext}"
                candidate = Path(chromhmm_dir) / pattern
                if candidate.exists():
                    bed_file = candidate
                    break
            
            if not bed_file:
                logger.warning(f"No chromHMM file found for EID {eid} in {chromhmm_dir}")
                continue

            #logger.info(f"Processing EID: {eid} | File: {bed_file}")
            logger.info(f"Processing EID: {eid} for {category_name} category")
            
            try:
                # Read first line to determine number of columns
                if str(bed_file).endswith('.gz'):
                    with gzip.open(bed_file, 'rt') as f:
                        # Skip the first line
                        f.readline()
                        # Read the second line to determine number of columns
                        second_line = f.readline().strip()
                        # logger.info(f"Second line of {bed_file}: {second_line}")
                        fields = second_line.split('\t')
                        num_cols = len(fields)
                        # logger.info(f"Number of columns in {bed_file}: {num_cols}")
                        
                        # Create column names
                        col_names = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + [f'extra_{i}' for i in range(num_cols - 6)]
                        # logger.info(f"Column names: {col_names}")
                        
                        # Read the file with correct number of columns, skipping the first line
                        df = pd.read_csv(f, sep='\t', header=None, names=col_names)
                else:
                    with open(bed_file, 'r') as f:
                        # Skip the first line
                        f.readline()
                        # Read the second line to determine number of columns
                        second_line = f.readline().strip()
                        # logger.info(f"Second line of {bed_file}: {second_line}")
                        fields = second_line.split('\t')
                        num_cols = len(fields)
                        # logger.info(f"Number of columns in {bed_file}: {num_cols}")
                        
                        # Create column names
                        col_names = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + [f'extra_{i}' for i in range(num_cols - 6)]
                        logger.info(f"Column names: {col_names}")
                        
                        # Read the file with correct number of columns, skipping the first line
                        df = pd.read_csv(bed_file, sep='\t', header=None, names=col_names, on_bad_lines='warn')
                
                logger.info(f"Successfully read {len(df)} rows from {bed_file}")
                # logger.info(f"Columns in dataframe: {df.columns.tolist()}")
                # logger.info(f"First few rows:\n{df.head()}")
                
                logger.info(f"Filtering for regulatory states: {states}")
                df = df[df['name'].isin(states)]
                logger.info(f"Found {len(df)} regions matching regulatory states")
                
                # Keep only the BED6 columns
                df = df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
                all_regions.append(df)
                
            except Exception as e:
                logger.error(f"Error processing file {bed_file}: {str(e)}")
                logger.error(f"Error type: {type(e)}")
                import traceback
                logger.error(f"Traceback: {traceback.format_exc()}")
                continue
        
        if not all_regions:
            logger.warning(f"No regions found for {category_name} category")
            return None
        
        # Combine all regions
        combined_regions = pd.concat(all_regions, ignore_index=True)
        
        # Save to BED file
        combined_regions.to_csv(output_file, sep='\t', header=False, index=False)
        
        print("")
        print("")
        logger.info(f"Total combined {len(combined_regions)} regions for {category_name} category")
        logger.info(f"Saved combined BED file to: {output_file}")
        print("")
        print("")
        return str(output_file)
        
    except Exception as e:
        logger.error(f"Error processing chromHMM files for {category_name}: {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        raise

def prepare_regulatory_regions(chromhmm_dir, output_dir, reg_states):
    """Prepare regulatory regions from chromHMM files.
    
    Args:
        chromhmm_dir (str): Directory containing chromHMM files
        output_dir (str): Directory to save processed files
        reg_states (dict): Dictionary containing regulatory states and EID categories
    
    Returns:
        dict: Dictionary mapping category names to their BED file paths
    """
    logger.info("Preparing regulatory regions from chromHMM files...")
    
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each category
        category_files = {}
        for category in ['fetal', 'stem_cell', 'adult']:
            bed_file = process_chromhmm_files(chromhmm_dir, category, output_dir, reg_states)
            category_files[category] = bed_file
            
        return category_files
        
    except Exception as e:
        logger.error(f"Error preparing regulatory regions: {str(e)}")
        raise


def save_results(results, output_dir):
    """Save analysis results to files."""
    logger.info("Saving results...")
    
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save summary statistics
        with open(output_dir / 'analysis_summary.txt', 'w') as f:
            f.write("cQTL Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total cQTLs: {results['total_regions']}\n")
            f.write(f"Overlapping cQTLs: {results['n_overlaps']}\n")
            f.write(f"Overlap fraction: {results['overlap_fraction']:.4f}\n\n")
            f.write("Fisher's Exact Test:\n")
            f.write(f"  Odds ratio: {results['fishers_exact']['odds_ratio']:.4f}\n")
            f.write(f"  P-value: {results['fishers_exact']['p_value']:.4e}\n\n")
            f.write("Permutation Test:\n")
            f.write(f"  Empirical P-value: {results['permutation_analysis']['empirical_p_value']:.4e}\n")
        
        # Save detailed results as JSON
        with open(output_dir / 'detailed_results.json', 'w') as f:
            json.dump(results, f, indent=2)
            
    except Exception as e:
        logger.error(f"Error saving results: {e}")
        raise



def sanity_check():
    """Perform sanity checks on main functions with provided file paths."""

    # Define test paths for cQTL analysis
    test_paths = {
        'base_dir': "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cQTL_H3K4me3_qvalBased_V2",
        'cqtl_file': "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/cfChIP_H3K4me3.combined.txt",
        'chromhmm_dir': "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/all.dense.browserFiles",
        'custom_random_bg': "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/WB.combined.txt",
        'genome_build': 'hg19',
        'regulatory_state': 'active',
        'mode': 'cqtl',
        'use_qvalue': True  # Add option to use q-value instead of p-value
    }

    # option for regulatory_states : ['bivalent', 'enhancer', 'promoter', 'active']

    try:
        print("\n=== Starting Sanity Checks ===\n")
        
        # Create base directory and regulatory state subdirectory
        base_dir = Path(test_paths['base_dir'])
        regulatory_dir = base_dir / test_paths["regulatory_state"]
        base_dir.mkdir(parents=True, exist_ok=True)
        regulatory_dir.mkdir(parents=True, exist_ok=True)
        
        # Update base_dir in test_paths to include regulatory state
        test_paths['base_dir'] = str(regulatory_dir)

        # 1. Test get_regulatory_states
        print("Testing get_regulatory_states()...")
        reg_states = get_regulatory_states(test_paths["regulatory_state"])

        # 2. Process input file based on mode
        print("\nProcessing input file...")
        if test_paths['mode'] == 'cqtl':
            print("Processing cQTL file...")
            # Use q-value if specified
            if test_paths.get('use_qvalue', False):
                input_bed = process_cqtl_file(
                    test_paths['cqtl_file'],
                    qval_threshold=0.05,  # Use q-value threshold
                    genome_build=test_paths['genome_build']
                )
            else:
                input_bed = process_cqtl_file(
                    test_paths['cqtl_file'],
                    pval_threshold=0.05,
                    genome_build=test_paths['genome_build']
                )
            # Convert DataFrame to BedTool
            input_bed_file = os.path.join(test_paths['base_dir'], "original_cqtls_with_id.bed")
            input_bed.to_csv(input_bed_file, sep="\t", index=False, header=False)
            input_bed = BedTool(input_bed_file)
            total_regions = len(input_bed)
        else:  # cfpeak mode
            print("Processing cfPeak file...")
            input_bed, total_regions = process_cfpeak_file(
                test_paths['cfpeak_file'],
                test_paths['base_dir']
            )
        print(f"✓ Processed {total_regions} regions")

        # 3. Test process_custom_random_bg
        print("\nTesting process_custom_random_bg()...")
        # Use q-value if specified
        if test_paths.get('use_qvalue', False):
            custom_bg, metadata, full_random_bed = process_custom_random_bg(
                test_paths['custom_random_bg'],
                test_paths['base_dir'],
                target_count=total_regions,
                qval_threshold=0.05,  # Use q-value threshold
                mode=test_paths['mode'],
                perform_permutation=True,
                n_permutations=100,
                return_full_bedtool=True
            )
        else:
            custom_bg, metadata, full_random_bed = process_custom_random_bg(
                test_paths['custom_random_bg'],
                test_paths['base_dir'],
                target_count=total_regions,
                pval_threshold=0.5,
                mode=test_paths['mode'],
                perform_permutation=True,
                n_permutations=100,
                return_full_bedtool=True
            )
        
        print(f"✓ Processed {len(custom_bg)} random background regions")
        if metadata['perform_permutation']:
            print(f"✓ Generated {metadata['n_permutations']} permutations")

        # 4. Test prepare_regulatory_regions
        print("\nTesting prepare_regulatory_regions()...")
        regulatory_regions = prepare_regulatory_regions(
            test_paths['chromhmm_dir'],
            test_paths['base_dir'],
            reg_states
        )
        regulatory_regions_dir = os.path.join(test_paths['base_dir'], "regulatory_regions")
        print(f"✓ Prepared regulatory regions for {len(regulatory_regions)} categories in regulatory dir: {regulatory_regions_dir}")

        # 5. Test define_specificity_categories
        print("\nTesting define_specificity_categories()...")
        category_files = define_specificity_categories(
            test_paths['base_dir'],
            regulatory_regions_dir,
            mode=test_paths['mode']  # Pass the mode parameter
        )
        print(f"✓ Created {len(category_files)} specificity categories")

        # 6. Test perform_overlap_analysis
        print("\nTesting perform_overlap_analysis()...")
        # Perform basic overlap analysis for each category
        for category_name, bed_file in category_files.items():
            print("")
            print(f"✓ Processing category: {category_name}")
            print("")
            custom_regions_bed = BedTool(bed_file)
            
            # Create source files dictionary
            source_files = {
                'input_file': test_paths['cqtl_file'] if test_paths['mode'] == 'cqtl' else test_paths['cfpeak_file'],
                'random_bg_file': test_paths['custom_random_bg'],
                'custom_regions_file': str(bed_file)
            }
            
            overlap_results = perform_overlap_analysis(
                input_bed,
                custom_regions_bed,
                (custom_bg, metadata, full_random_bed),
                test_paths['base_dir'],
                category_name,
                mode=test_paths['mode'],
                source_files=source_files,
                perform_permutation=True,
                n_permutations=100
            )

        print(f"✓ Found {overlap_results['n_overlaps']} overlaps")
        print(f"✓ Overlap fraction: {overlap_results['overlap_fraction']:.4f}")

    except Exception as e:
        print(f"\n❌ Error during sanity check: {str(e)}")
        raise


if __name__ == "__main__":

    # Run sanity checks
    sanity_check()



#-----------------------------------------------------------------------------------------

# One time run only

#-----------------------------------------------------------------------------------------


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


#-----------------------------------------------------------------------------------------

# One time run only

#-----------------------------------------------------------------------------------------


print("Processing cQTL file...")

# Define test paths for cQTL analysis
test_paths = {
    'base_dir': "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2",
    'cqtl_file': "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/cfChIP_H3K4me3.combined.txt",
    'custom_random_bg': "/Users/chhetribsurya/sc1238/datasets/projects/cfChIP_project/cancer_developmental_region_analysis/WB.combined.txt",
    'genome_build': 'hg19',
}

# Create base directory if it doesn't exist
base_dir = test_paths['base_dir']
if not os.path.exists(base_dir):
    os.makedirs(base_dir, exist_ok=True)
    print(f"✓ Created base directory: {base_dir}")
else:
    print(f"✓ Base directory already exists: {base_dir}")


input_bed = process_cqtl_file(
        test_paths['cqtl_file'],
        qval_threshold=0.05,  # Use q-value threshold
        genome_build=test_paths['genome_build']
    )

input_bed_random = process_cqtl_file(
        test_paths['custom_random_bg'],
        qval_threshold=0.05,  # Use q-value threshold
        genome_build=test_paths['genome_build']
    )

# Convert DataFrame to BedTool
input_bed_file = os.path.join(test_paths['base_dir'], "original_cqtls_with_id.bed")
input_bed.to_csv(input_bed_file, sep="\t", index=False, header=False)

input_bed_random_file = os.path.join(test_paths['base_dir'], "random_cqtls_with_id.bed")
input_bed_random.to_csv(input_bed_random_file, sep="\t", index=False, header=False)

# Create BedTool objects
input_bedtool = BedTool(input_bed_file)
input_bed_randomtool = BedTool(input_bed_random_file)

# Find input regions NOT overlapping with random background
input_unique_not_in_random = input_bedtool.intersect(input_bed_randomtool, v=True)
unique_outfile = os.path.join(test_paths['base_dir'], "cfcQTLs_unique_to_WBCcQTLs.bed")
input_unique_not_in_random.saveas(unique_outfile)
print(f"Saved unique input cQTLs not in random background to: {unique_outfile}")

# Find input regions NOT overlapping with random background
input_unique_not_in_random = input_bed_randomtool.intersect(input_bedtool, v=True)
unique_outfile = os.path.join(test_paths['base_dir'], "WBCcQTLs_unique_to_cfcQTLs.bed")
input_unique_not_in_random.saveas(unique_outfile)
print(f"Saved unique input cQTLs not in random background to: {unique_outfile}")

#-----------------------------------------------------------------------------------------
# --- CONFIGURATION ---

input_bed = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2/WBCcQTLs_unique_to_cfcQTLs.bed"
input_bed = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2/cfcQTLs_unique_to_WBCcQTLs.bed"

input_bed = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2/random_cqtls_with_id.bed"
input_bed = "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2/original_cqtls_with_id.bed"

# Read the input bed file
input_bed = BedTool(input_bed)

base_dir = test_paths['base_dir']
base_dir="/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/for_david/NewTest_cQTL_H3K4me3_qvalBased_V2"

cat_dir = os.path.join(base_dir, "active/regulatory_regions/specificity_categories")

out_dir = os.path.join(base_dir, "WBCcQTLs_unique_to_cfcQTLs")
out_dir = os.path.join(base_dir, "cfcQTLs_unique_to_WBCcQTLs")

out_dir = os.path.join(base_dir, "full_original_WBCcQTLs")
out_dir = os.path.join(base_dir, "full_original_cfcQTLs")

# --- CONFIGURATION ---

os.makedirs(out_dir, exist_ok=True)
analysis_dir = os.path.join(out_dir, "analysis")
os.makedirs(analysis_dir, exist_ok=True)

# --- LOAD CATEGORY FILES ---
cat_files = [f for f in os.listdir(cat_dir) if f.endswith('.bed')]
cat_files.sort()  # for consistent order
cat_labels = [os.path.splitext(f)[0].replace("_", " ") for f in cat_files]  # Clean for plot
cat_paths = [os.path.join(cat_dir, f) for f in cat_files]

# --- INTERSECT AND SAVE ---
intersected = {}
for fname, path in zip(cat_files, cat_paths):
    bed = BedTool(path)
    intersect = input_bed.intersect(bed, u=True)
    out_path = os.path.join(out_dir, f"{fname.replace('.bed','')}_intersected.bed")
    intersect.saveas(out_path)
    intersected[fname] = set(f"{x.chrom}:{x.start}-{x.end}" for x in intersect)
    print(f"Saved {len(intersected[fname])} unique overlaps for {fname} to {out_path}")

# --- BUILD MEMBERSHIP TABLE FOR UPSET ---
all_ids = set.union(*intersected.values())
membership = []
for cid in all_ids:
    present = [fname for fname, ids in intersected.items() if cid in ids]
    membership.append(present)
upset_data = from_memberships(membership)
# Aggregate to ensure unique index
upset_data = upset_data.groupby(level=list(range(upset_data.index.nlevels))).sum()

# Remove '.bed' from index names for prettier labels
upset_data.index.names = [name.replace('.bed', '').replace('_', ' ').title() for name in upset_data.index.names]

# --- SAVE SUMMARY TABLE ---
summary_df = upset_data.reset_index()
summary_df.columns = cat_files + ['count']
summary_df.to_csv(os.path.join(analysis_dir, "upset_summary.csv"), index=False)

# --- PLOT UPSET ---
import matplotlib.pyplot as plt
from matplotlib import cm

plt.figure(figsize=(10, 7))
up = UpSet(
    upset_data,
    show_counts=True,
    sort_by='cardinality',
    element_size=None
)
up.plot()

# Get the intersection size axes robustly
ax = plt.gca()
ax.grid(False)

# Use a colormap for freestyle coloring
cmap = cm.get_cmap('tab10')
num_bars = len(ax.patches)
colors = [cmap(i % cmap.N) for i in range(num_bars)]

# Color each bar with a different color
for patch, color in zip(ax.patches, colors):
    patch.set_facecolor(color)

# Reduce font size for all labels
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(10)
ax.title.set_fontsize(11)
ax.xaxis.label.set_fontsize(11)
ax.yaxis.label.set_fontsize(11)

plt.suptitle("UpSet Plot of WBC-cQTLs Overlapping Specificity Categories", fontsize=13)
plt.suptitle("UpSet Plot of cfcQTLs Overlapping Specificity Categories", fontsize=13)
plt.savefig(os.path.join(analysis_dir, "upset_plot.pdf"), bbox_inches='tight')
plt.close()

# --- PLOT VENN (if <=3 categories) ---
if len(cat_files) == 2:
    sets = [intersected[cat_files[0]], intersected[cat_files[1]]]
    plt.figure(figsize=(6, 6))
    venn2(sets, set_labels=[cat_files[0].replace("_", "").replace(".bed", ""), cat_files[1].replace("_", "").replace(".bed", "")])
    plt.title("Venn Diagram of cfcQTLs Overlap")
    plt.savefig(os.path.join(analysis_dir, "venn2_plot.pdf"), bbox_inches='tight')
    plt.close()
elif len(cat_files) == 3:
    sets = [intersected[cat_files[0]], intersected[cat_files[1]], intersected[cat_files[2]]]
    plt.figure(figsize=(7, 7))
    venn3(sets, set_labels=[cat_files[0].replace("_", "").replace(".bed", ""), cat_files[1].replace("_", "").replace(".bed", ""), cat_files[2].replace("_", "").replace(".bed", "")])
    plt.title("Venn Diagram of cfcQTLs Overlap")
    plt.savefig(os.path.join(analysis_dir, "venn3_plot.pdf"), bbox_inches='tight')
    plt.close()

print(f"All intersections and plots saved in {analysis_dir}")


#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------

# --- BARPLOT OF UNIQUE OVERLAPS FOR EACH CATEGORY (FIXED ORDER & COLORS, CLEAN Y-TICKS) ---

# Fixed color mapping for each category
category_color_map = {
    "fetal_specific.bed": "#F1C40F",
    "developmental_specific.bed": "#E74C3C",
    "stem_cell_specific.bed": "#3498DB",
    "adult_specific.bed": "#2ECC71"
}

# Prepare data in the fixed order
fixed_order = [
    "fetal_specific.bed",
    "developmental_specific.bed",
    "stem_cell_specific.bed",
    "adult_specific.bed"
]
bar_names = [name.replace('.bed', '').replace('_', ' ').title() for name in fixed_order]
bar_counts = [len(intersected[name]) for name in fixed_order]
bar_colors = [category_color_map[name] for name in fixed_order]

fig, ax = plt.subplots(figsize=(7, 6))
bars = ax.bar(range(len(bar_counts)), bar_counts,
              color=bar_colors,
              width=0.6, edgecolor='black')

# X labels (smaller font size)
ax.set_xticks(range(len(bar_counts)))
ax.set_xticklabels(bar_names, rotation=20, ha='right')
ax.tick_params(axis='x', labelsize=9)  # This will always work

# Legends
legend_elements = [Rectangle((0, 0), 1, 1, facecolor=category_color_map[name], edgecolor='black', label=bar_names[i])
                   for i, name in enumerate(fixed_order)]
ax.legend(
    handles=legend_elements,
    loc='upper left',
    frameon=False,
    handleheight=1.5,
    handlelength=1.5,
    borderpad=0.8,
    labelspacing=0.8,
    fontsize=10
)

# Despine and offset
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
ax.spines['left'].set_position(('outward', 8))
ax.spines['bottom'].set_position(('outward', 8))
ax.tick_params(axis='both', which='both', length=6, width=1, direction='out')
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_ylabel('Number of Unique Overlaps', fontsize=13)

# Y-axis ticks (fix duplication/overlap and round up)
y_min = 0
raw_y_max = max(bar_counts)
# For large numbers, round up to the next 100 or 1000 for a clean axis
if raw_y_max > 1000:
    y_max = int(np.ceil(raw_y_max / 100.0)) * 100
    step = 1000 if y_max > 5000 else 500
elif raw_y_max > 100:
    y_max = int(np.ceil(raw_y_max / 10.0)) * 10
    step = 100
else:
    y_max = int(np.ceil(raw_y_max))
    step = 1

yticks = np.arange(y_min, y_max + step, step)
ax.set_ylim([y_min, y_max])
ax.set_yticks(yticks)
ax.set_yticklabels([str(int(y)) for y in yticks])

# Add value labels
for i, count in enumerate(bar_counts):
    ax.text(i, count + y_max*0.02, f"{count}", ha='center', va='bottom', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, 'unique_overlap_counts_barplot.pdf'), bbox_inches='tight')
plt.close()





