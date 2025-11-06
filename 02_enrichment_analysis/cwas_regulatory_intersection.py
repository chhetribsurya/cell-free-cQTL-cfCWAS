#!/usr/bin/env python3
"""
cfCWAS Regulatory Element Intersection Analysis Script (pybedtools version)

This script analyzes the intersection between cfCWAS (cell-free Cistrome-Wide Association Study) 
loci and various regulatory element files (adult, developmental, fetal, and stem cell specific) 
to understand the relationship between pleiotropy, trait associations, and regulatory context.

Uses pybedtools for genomic intersections and provides step-by-step execution for debugging.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Import pybedtools
try:
    import pybedtools
    print("pybedtools imported successfully")
except ImportError:
    print("Error: pybedtools not found. Please install with: pip install pybedtools")
    exit(1)

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def convert_to_bed(id_str):
    """Convert cfCWAS ID format to BED format"""
    # Format: "chr10:32343900-32345900" -> "chr10\t32343900\t32345900"
    match = re.match(r'(chr[^:]+):(\d+)-(\d+)', id_str)
    if match:
        chrom, start, end = match.groups()
        return f"{chrom}\t{start}\t{end}"
    return None

def load_cwas_data(cwas_file, output_dir):
    """Load and preprocess cfCWAS data"""
    print("STEP 1: Loading cfCWAS data...")
    
    try:
        # Read cfCWAS file
        cwas_data = pd.read_csv(cwas_file)
        print(f"Successfully loaded {len(cwas_data)} cfCWAS loci from {cwas_file}")
        
        # Convert ID column to BED format for intersection
        cwas_data['bed_coords'] = cwas_data['ID'].apply(convert_to_bed)
        
        # Check for any failed conversions
        failed_conversions = cwas_data['bed_coords'].isna().sum()
        if failed_conversions > 0:
            print(f"Warning: {failed_conversions} IDs could not be converted to BED format")
        
        # Extract chromosome and position information
        cwas_data[['chrom', 'start', 'end']] = cwas_data['bed_coords'].str.split('\t', expand=True)
        cwas_data['start'] = cwas_data['start'].astype(int)
        cwas_data['end'] = cwas_data['end'].astype(int)
        
        # Create temporary BED file for intersection
        cwas_bed_file = f"{output_dir}/cfcwas_loci.bed"
        cwas_data[['chrom', 'start', 'end', 'ID']].to_csv(
            cwas_bed_file, sep='\t', header=False, index=False
        )
        
        print(f"Created BED file: {cwas_bed_file}")
        print(f"Pleiotropy distribution:\n{cwas_data['pleiotropy'].value_counts()}")
        
        return cwas_data, cwas_bed_file
        
    except Exception as e:
        print(f"Error loading cfCWAS data: {e}")
        raise

def run_pybedtools_intersection(file_a, file_b, output_file):
    """Run pybedtools intersection with -wa -wb flags"""
    print(f"  Running intersection: {os.path.basename(file_a)} ∩ {os.path.basename(file_b)}")
    
    try:
        # Create BedTool objects
        bed_a = pybedtools.BedTool(file_a)
        bed_b = pybedtools.BedTool(file_b)
        
        # Perform intersection
        intersection = bed_a.intersect(bed_b, wa=True, wb=True)
        
        # Save to file
        intersection.saveas(output_file)
        
        # Get intersection count
        intersection_count = len(intersection)
        print(f"  Found {intersection_count} intersections")
        
        return True, intersection_count
        
    except Exception as e:
        print(f"  Error in pybedtools intersection: {e}")
        return False, 0

def analyze_intersection(reg_type, intersection_file, cwas_data, output_dir):
    """Analyze intersection results for a specific regulatory type"""
    print(f"  Analyzing {reg_type} intersection results...")
    
    try:
        if not os.path.exists(intersection_file) or os.path.getsize(intersection_file) == 0:
            print(f"  No intersections found for {reg_type}")
            intersection_results = pd.DataFrame()
            intersection_stats = {
                'total_intersections': 0,
                'unique_cwas_loci': 0,
                'unique_cQTL_elements': 0
            }
            return intersection_results, intersection_stats
        
        # Read intersection file using pybedtools
        intersection_bed = pybedtools.BedTool(intersection_file)
        
        # Convert to DataFrame
        intersection_data = []
        for interval in intersection_bed:
            intersection_data.append({
                'cwas_chrom': interval.chrom,
                'cwas_start': interval.start,
                'cwas_end': interval.end,
                'cwas_id': interval.name,
                'cQTL_chrom': interval.fields[4],
                'cQTL_start': int(interval.fields[5]),
                'cQTL_end': int(interval.fields[6]),
                'cQTL_id': interval.fields[7]
            })
        
        intersection_df = pd.DataFrame(intersection_data)
        
        # Merge with original CWAS data
        merged_data = intersection_df.merge(
            cwas_data, left_on='cwas_id', right_on='ID', how='left'
        )
        
        # Calculate statistics
        unique_cwas = merged_data['cwas_id'].nunique()
        unique_cQTL = merged_data['cQTL_id'].nunique()
        total_intersections = len(merged_data)
        
        intersection_stats = {
            'total_intersections': total_intersections,
            'unique_cwas_loci': unique_cwas,
            'unique_cQTL_elements': unique_cQTL
        }
        
        print(f"  {reg_type}: {total_intersections} intersections, "
              f"{unique_cwas} unique cfCWAS loci, {unique_cQTL} unique cQTL elements")
        
        # Save detailed intersection results
        merged_data.to_csv(
            f"{output_dir}/intersections/cfcwas_{reg_type}_detailed.csv", index=False
        )
        print(f"  Saved detailed results to: {output_dir}/intersections/cfcwas_{reg_type}_detailed.csv")
        
        return merged_data, intersection_stats
        
    except Exception as e:
        print(f"  Error analyzing {reg_type} intersection: {e}")
        raise

def perform_intersections(cwas_bed_file, regulatory_dir, output_dir, cwas_data):
    """Perform intersections with all regulatory element files"""
    print("\nSTEP 2: Performing intersections with regulatory elements...")
    
    # Regulatory element file patterns
    regulatory_files = {
        'adult_specific': 'adult_specific_intersected.bed',
        'developmental_specific': 'developmental_specific_intersected.bed', 
        'fetal_specific': 'fetal_specific_intersected.bed',
        'stem_cell_specific': 'stem_cell_specific_intersected.bed'
    }
    
    intersection_results = {}
    intersection_stats = {}
    
    for reg_type, bed_file in regulatory_files.items():
        reg_file_path = os.path.join(regulatory_dir, bed_file)
        
        if not os.path.exists(reg_file_path):
            print(f"Warning: {reg_file_path} not found, skipping {reg_type}...")
            continue
            
        print(f"\nProcessing {reg_type} regulatory elements...")
        
        # Perform intersection using pybedtools
        intersection_file = f"{output_dir}/intersections/cfcwas_{reg_type}_intersection.bed"
        success, count = run_pybedtools_intersection(cwas_bed_file, reg_file_path, intersection_file)
        
        if success:
            # Analyze intersection results
            results, stats = analyze_intersection(reg_type, intersection_file, cwas_data, output_dir)
            intersection_results[reg_type] = results
            intersection_stats[reg_type] = stats
        else:
            print(f"Failed to process {reg_type} regulatory elements")
    
    print(f"\nCompleted intersections for {len(intersection_stats)} regulatory types")
    return intersection_results, intersection_stats

def clean_trait_name(trait):
    """Clean trait name for display by removing common prefixes and formatting"""
    clean_trait = trait.replace('PASS_', '').replace('_', ' ')
    # Remove "Non-cancer illness code, self-reported: " prefix
    if clean_trait.startswith('Non-cancer illness code, self-reported: '):
        clean_trait = clean_trait.replace('Non-cancer illness code, self-reported: ', '')
    return clean_trait

def is_non_ukb_trait(trait):
    """Check if a trait is non-ukb related"""
    trait_lower = trait.lower()
    non_ukb_keywords = ['non-ukb', 'non_ukb', 'non ukb', 'nonukb']
    return any(keyword in trait_lower for keyword in non_ukb_keywords)

def analyze_traits(intersected_data):
    """Analyze trait associations in intersected data"""
    trait_counter = Counter()
    
    # Check if all_descriptions column exists
    if 'all_descriptions' in intersected_data.columns:
    for descriptions in intersected_data['all_descriptions'].dropna():
            if descriptions and descriptions != '':
        traits = descriptions.split('; ')
        trait_counter.update(traits)
    
    return trait_counter.most_common()

def prioritize_traits_for_plotting(trait_counts, max_traits=35):
    """
    Prioritize traits for plotting with special focus on non-ukb and cancer-related traits
    
    Args:
        trait_counts: List of (trait, count) tuples from Counter.most_common()
        max_traits: Maximum number of traits to include in plot
    
    Returns:
        List of prioritized (trait, count) tuples
    """
    if len(trait_counts) <= max_traits:
        return trait_counts
    
    # Separate traits into priority categories
    cancer_traits = []
    non_ukb_traits = []
    other_traits = []
    
    for trait, count in trait_counts:
        trait_lower = trait.lower()
        
        # Check for cancer-related terms
        cancer_keywords = ['cancer', 'carcinoma', 'tumor', 'neoplasm', 'malignant', 'oncology', 
                          'lung cancer', 'breast cancer', 'prostate cancer', 'colorectal cancer',
                          'ovarian cancer', 'pancreatic cancer', 'liver cancer', 'brain cancer',
                          'leukemia', 'lymphoma', 'melanoma', 'sarcoma']
        
        # Check for non-ukb terms
        non_ukb_keywords = ['non-ukb', 'non_ukb', 'non ukb', 'nonukb']
        
        if any(keyword in trait_lower for keyword in cancer_keywords):
            cancer_traits.append((trait, count))
        elif any(keyword in trait_lower for keyword in non_ukb_keywords):
            non_ukb_traits.append((trait, count))
        else:
            other_traits.append((trait, count))
    
    # Combine prioritized traits
    prioritized_traits = []
    
    # First add all cancer-related traits
    prioritized_traits.extend(cancer_traits)
    
    # Then add all non-ukb traits
    prioritized_traits.extend(non_ukb_traits)
    
    # Fill remaining slots with other traits (sorted by count)
    remaining_slots = max_traits - len(prioritized_traits)
    if remaining_slots > 0:
        prioritized_traits.extend(other_traits[:remaining_slots])
    
    # If we still have more than max_traits, take the top ones
    if len(prioritized_traits) > max_traits:
        prioritized_traits = prioritized_traits[:max_traits]
    
    return prioritized_traits

def generate_statistics(intersection_stats, cwas_data, intersection_results, output_dir):
    """Generate comprehensive statistics"""
    print("\nSTEP 3: Generating comprehensive statistics...")
    
    try:
        stats_data = []
        
        for reg_type, stats in intersection_stats.items():
            # Basic intersection stats
            row = {
                'regulatory_type': reg_type,
                'total_intersections': stats['total_intersections'],
                'unique_cwas_loci': stats['unique_cwas_loci'],
                'unique_cQTL_elements': stats['unique_cQTL_elements'],
                'total_cwas_loci': len(cwas_data),
                'intersection_rate': stats['unique_cwas_loci'] / len(cwas_data) * 100
            }
            
            # Pleiotropy analysis if intersections exist
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                intersected_data = intersection_results[reg_type]
                
                # Pleiotropy distribution
                pleiotropy_counts = intersected_data['pleiotropy'].value_counts()
                row.update({
                    'trait_specific_count': pleiotropy_counts.get('trait-specific', 0),
                    'pleiotropic_count': pleiotropy_counts.get('pleiotropic', 0),
                    'trait_specific_rate': pleiotropy_counts.get('trait-specific', 0) / stats['unique_cwas_loci'] * 100,
                    'pleiotropic_rate': pleiotropy_counts.get('pleiotropic', 0) / stats['unique_cwas_loci'] * 100
                })
                
                # Trait analysis
                trait_counts = analyze_traits(intersected_data)
                row['top_traits'] = '; '.join([f"{trait}({count})" for trait, count in trait_counts[:5]])
            
            stats_data.append(row)
        
        # Create statistics DataFrame
        stats_df = pd.DataFrame(stats_data)
        stats_df.to_csv(f"{output_dir}/statistics/intersection_statistics.csv", index=False)
        
        print(f"Saved statistics to: {output_dir}/statistics/intersection_statistics.csv")
        
        # Print summary
        print("\nIntersection Statistics Summary:")
        print(stats_df.to_string(index=False))
        
        return stats_df
        
    except Exception as e:
        print(f"Error generating statistics: {e}")
        raise

def plot_intersection_overview(intersection_stats, cwas_data, output_dir):
    """Plot intersection overview"""
    print("  Creating intersection overview plot...")
    
    try:
        # Set publication-ready style
        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), facecolor='white')
        
        # Intersection counts
        reg_types = list(intersection_stats.keys())
        intersection_counts = [intersection_stats[rt]['total_intersections'] for rt in reg_types]
        
        # Consistent color mapping for each regulatory type - all dark pink
        color_map = {
            'adult_specific': '#E91E63',      # Dark pink
            'developmental_specific': '#E91E63',  # Dark pink
            'fetal_specific': '#E91E63',      # Dark pink
            'stem_cell_specific': '#E91E63'   # Dark pink
        }
        
        # Get colors for each regulatory type
        colors = [color_map.get(rt, '#E91E63') for rt in reg_types]
        
        x = np.arange(len(reg_types))
        
        # Plot 1: Intersection counts only
        ax1.bar(x, intersection_counts, width=0.6, alpha=0.8, color=colors, 
                       edgecolor='black', linewidth=0.5)
        
        ax1.set_xlabel('Specific Regulatory Element Type', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Peak Count', fontsize=12, fontweight='bold')
        ax1.set_title('cfCWAS-Specific Regulatory Element Intersections', fontsize=14, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels([f"{rt.replace('_', ' ').title()}" for rt in reg_types], rotation=45, ha='right')
        ax1.grid(True, alpha=0.3, linestyle='--')
        
        # Set y-axis limits and ticks
        max_count = max(intersection_counts)
        ax1.set_ylim(0, max_count * 1.1)
        ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        
        # Despine with offset
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_linewidth(1.5)
        ax1.spines['bottom'].set_linewidth(1.5)
        ax1.spines['left'].set_position(('outward', 10))
        ax1.spines['bottom'].set_position(('outward', 10))
        
        # Intersection rates
        total_cwas = len(cwas_data)
        intersection_rates = [intersection_stats[rt]['unique_cwas_loci'] / total_cwas * 100 
                            for rt in reg_types]
        
        # Plot 2: Intersection rates
        ax2.bar(reg_types, intersection_rates, alpha=0.8, color=colors, 
                       edgecolor='black', linewidth=0.5)
        ax2.set_xlabel('Specific Regulatory Element Type', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Intersection Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_title('cfCWAS Loci Intersection Rate by Specific Regulatory Type', fontsize=14, fontweight='bold')
        ax2.set_xticks(range(len(reg_types)))
        ax2.set_xticklabels([f"{rt.replace('_', ' ').title()}" for rt in reg_types], rotation=45, ha='right')
        ax2.grid(True, alpha=0.3, linestyle='--')
        
        # Set y-axis limits and ticks (allow decimals for percentages)
        max_rate = max(intersection_rates)
        ax2.set_ylim(0, max_rate * 1.1)
        # Allow decimals for percentage values
        ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2f}'))
        
        # Despine with offset
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_linewidth(1.5)
        ax2.spines['bottom'].set_linewidth(1.5)
        ax2.spines['left'].set_position(('outward', 10))
        ax2.spines['bottom'].set_position(('outward', 10))
        
        plt.tight_layout()
        fig.canvas.draw()  # Ensure ticks are up to date

        # For ax1 (left panel)
        x_ticks = ax1.get_xticks()
        y_ticks = ax1.get_yticks()
        if len(x_ticks) > 1:
            ax1.set_xlim(x_ticks[0], x_ticks[-1])
        if len(y_ticks) > 1:
            ax1.set_ylim(y_ticks[0], y_ticks[-1])
        fig.canvas.draw()
        if len(x_ticks) > 1:
            ax1.spines['bottom'].set_bounds(x_ticks[0], x_ticks[-1])
        if len(y_ticks) > 1:
            ax1.spines['left'].set_bounds(y_ticks[0], y_ticks[-1])

        # For ax2 (right panel)
        x_ticks = ax2.get_xticks()
        y_ticks = ax2.get_yticks()
        if len(x_ticks) > 1:
            ax2.set_xlim(x_ticks[0], x_ticks[-1])
        if len(y_ticks) > 1:
            ax2.set_ylim(y_ticks[0], y_ticks[-1])
        fig.canvas.draw()
        if len(x_ticks) > 1:
            ax2.spines['bottom'].set_bounds(x_ticks[0], x_ticks[-1])
        if len(y_ticks) > 1:
            ax2.spines['left'].set_bounds(y_ticks[0], y_ticks[-1])
        
        plt.savefig(f"{output_dir}/plots/intersection_overview.png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(f"{output_dir}/plots/intersection_overview.pdf", bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        print("  Saved intersection overview plot")
        
    except Exception as e:
        print(f"  Error creating intersection overview plot: {e}")
        raise

def plot_pleiotropy_analysis(intersection_results, output_dir):
    """Plot pleiotropy analysis"""
    print("  Creating pleiotropy analysis plots...")
    
    try:
        # Set publication-ready style
        plt.style.use('default')
        fig, axes = plt.subplots(2, 2, figsize=(15, 12), facecolor='white')
        axes = axes.ravel()
        
        regulatory_files = ['adult_specific', 'developmental_specific', 'fetal_specific', 'stem_cell_specific']
        
        # Define color scheme for legend
        color_scheme = {
            'trait-specific': '#DC143C',  # Crimson red
            'pleiotropic': '#4169E1',     # Royal blue
            'other': '#808080'            # Gray for any other types
        }
        
        for i, reg_type in enumerate(regulatory_files):
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                data = intersection_results[reg_type]
                pleiotropy_counts = data['pleiotropy'].value_counts()
                
                # Define elegant color scheme: red for trait-specific, blue for pleiotropic
                pie_colors = []
                for pleiotropy_type in pleiotropy_counts.index:
                    if pleiotropy_type == 'trait-specific':
                        pie_colors.append(color_scheme['trait-specific'])
                    elif pleiotropy_type == 'pleiotropic':
                        pie_colors.append(color_scheme['pleiotropic'])
                    else:
                        pie_colors.append(color_scheme['other'])
                
                # Pie chart with publication-ready styling
                wedges, texts, autotexts = axes[i].pie(pleiotropy_counts.values, 
                                                      labels=pleiotropy_counts.index, 
                                                      autopct='%1.1f%%',
                                                      colors=pie_colors,
                                                      startangle=90,
                                                      textprops={'fontsize': 10, 'fontweight': 'bold'})
                
                # Style the percentage text for better contrast
                for autotext in autotexts:
                    autotext.set_color('black')
                    autotext.set_fontweight('bold')
                
                axes[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements\nPleiotropy Distribution', 
                                fontsize=12, fontweight='bold', pad=20)
            else:
                axes[i].text(0.5, 0.5, 'No intersections', ha='center', va='center', 
                           transform=axes[i].transAxes, fontsize=14, fontweight='bold')
                axes[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements', 
                                fontsize=12, fontweight='bold')
            
            # Remove spines for pie charts
            axes[i].spines['top'].set_visible(False)
            axes[i].spines['right'].set_visible(False)
            axes[i].spines['left'].set_visible(False)
            axes[i].spines['bottom'].set_visible(False)
        
        # Create a single legend for all pie charts
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor=color_scheme['trait-specific'], 
                         edgecolor='black', linewidth=0.5, label='Trait-specific'),
            plt.Rectangle((0, 0), 1, 1, facecolor=color_scheme['pleiotropic'], 
                         edgecolor='black', linewidth=0.5, label='Pleiotropic'),
            plt.Rectangle((0, 0), 1, 1, facecolor=color_scheme['other'], 
                         edgecolor='black', linewidth=0.5, label='Other')
        ]
        
        # Add legend to the figure (positioned at the bottom center)
        fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.05), 
                  ncol=3, fontsize=12, frameon=True, fancybox=True, shadow=True)
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)  # Make room for the legend
        plt.savefig(f"{output_dir}/plots/pleiotropy_analysis.png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(f"{output_dir}/plots/pleiotropy_analysis.pdf", bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        print("  Saved pleiotropy analysis plots")
        
    except Exception as e:
        print(f"  Error creating pleiotropy analysis plots: {e}")
        raise

def plot_trait_analysis(intersection_results, output_dir):
    """Plot trait association analysis"""
    print("  Creating trait analysis plots...")
    
    try:
        # Set publication-ready style
        plt.style.use('default')
        
        # Create two separate figures: one for top 10 traits and one for all traits
        regulatory_files = ['adult_specific', 'developmental_specific', 'fetal_specific', 'stem_cell_specific']
        # Consistent color mapping for each regulatory type - all dark pink
        color_map = {
            'adult_specific': '#E91E63',      # Dark pink
            'developmental_specific': '#E91E63',  # Dark pink
            'fetal_specific': '#E91E63',      # Dark pink
            'stem_cell_specific': '#E91E63'   # Dark pink
        }
        
        # Create individual plots for each regulatory region
        for reg_type in regulatory_files:
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                data = intersection_results[reg_type]
                trait_counts = analyze_traits(data)
                
                if trait_counts:
                    # Create individual plot for this regulatory type
                    create_individual_trait_plot(reg_type, trait_counts, color_map[reg_type], output_dir)
        
        # Create combined plots (2x2 subplot)
        create_combined_trait_plots(intersection_results, color_map, output_dir)
        
        print("  Saved trait analysis plots (individual and combined)")
        
    except Exception as e:
        print(f"  Error creating trait analysis plots: {e}")
        raise

def create_individual_trait_plot(reg_type, trait_counts, color, output_dir):
    """Create individual trait plot for a specific regulatory type"""
    try:
        # Use prioritized traits (max 35)
        prioritized_traits = prioritize_traits_for_plotting(trait_counts, max_traits=35)
        
        # Sort by count (descending order)
        prioritized_traits = sorted(prioritized_traits, key=lambda x: x[1], reverse=True)
        
        traits, counts = zip(*prioritized_traits)
        
        # Clean trait names for better display
        clean_traits = [clean_trait_name(trait) for trait in traits]
        
        # Ensure counts are integers
        counts = [int(count) for count in counts]
        
        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(12, max(8, len(traits) * 0.3)), facecolor='white')
        
        # Create horizontal bar plot
        y_positions = np.arange(len(clean_traits))
        ax.barh(y_positions, counts, height=0.6, alpha=0.8, 
                color=color, edgecolor='black', linewidth=0.5)
        ax.set_yticks(y_positions)
        
        # Create y-axis labels with highlighting for non-ukb traits
        y_labels = []
        for i, (trait, clean_trait) in enumerate(zip(traits, clean_traits)):
            if is_non_ukb_trait(trait):
                y_labels.append(f"★ {clean_trait}")
            else:
                y_labels.append(clean_trait)
        
        ax.set_yticklabels(y_labels, fontsize=10)
        
        # Color non-ukb trait labels red
        for i, (trait, label) in enumerate(zip(traits, y_labels)):
            if is_non_ukb_trait(trait):
                ax.get_yticklabels()[i].set_color('red')
                ax.get_yticklabels()[i].set_fontweight('bold')
        ax.set_xlabel('Peak Count', fontsize=12, fontweight='bold')
        ax.set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements\nTop {len(traits)} Trait Associations', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        # Set x-axis limits and ticks
        max_count = max(counts)
        ax.set_xlim(0, max_count * 1.1)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        
        # Despine with offset
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_position(('outward', 10))
        ax.spines['bottom'].set_position(('outward', 10))
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/trait_analysis_{reg_type}_individual.png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(f"{output_dir}/plots/trait_analysis_{reg_type}_individual.pdf", bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        print(f"    Saved individual plot for {reg_type}")
        
    except Exception as e:
        print(f"    Error creating individual plot for {reg_type}: {e}")

def create_combined_trait_plots(intersection_results, color_map, output_dir):
    """Create combined trait plots (2x2 subplot)"""
    try:
        regulatory_files = ['adult_specific', 'developmental_specific', 'fetal_specific', 'stem_cell_specific']
        
        # Figure 1: Top 10 traits (2x2 subplot)
        fig1, axes1 = plt.subplots(2, 2, figsize=(20, 12), facecolor='white')
        axes1 = axes1.ravel()
        
        for i, reg_type in enumerate(regulatory_files):
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                data = intersection_results[reg_type]
                trait_counts = analyze_traits(data)
                
                if trait_counts:
                    # Use prioritized traits (max 10 for combined plot)
                    prioritized_traits = prioritize_traits_for_plotting(trait_counts, max_traits=10)
                    
                    # Sort by count (descending order)
                    prioritized_traits = sorted(prioritized_traits, key=lambda x: x[1], reverse=True)
                    
                    traits, counts = zip(*prioritized_traits)
                    
                    # Clean trait names for better display
                    clean_traits = [clean_trait_name(trait) for trait in traits]
                    
                    # Ensure counts are integers
                    counts = [int(count) for count in counts]
                    
                    # Create horizontal bar plot with consistent bar width
                    y_positions = np.arange(len(clean_traits))
                    axes1[i].barh(y_positions, counts, height=0.6, alpha=0.8, 
                                 color=color_map[reg_type], edgecolor='black', linewidth=0.5)
                    axes1[i].set_yticks(y_positions)
                    
                    # Create y-axis labels with highlighting for non-ukb traits
                    y_labels = []
                    for j, (trait, clean_trait) in enumerate(zip(traits, clean_traits)):
                        if is_non_ukb_trait(trait):
                            y_labels.append(f"★ {clean_trait}")
                        else:
                            y_labels.append(clean_trait)
                    
                    axes1[i].set_yticklabels(y_labels, fontsize=10)
                    
                    # Color non-ukb trait labels red
                    for j, (trait, label) in enumerate(zip(traits, y_labels)):
                        if is_non_ukb_trait(trait):
                            axes1[i].get_yticklabels()[j].set_color('red')
                            axes1[i].get_yticklabels()[j].set_fontweight('bold')
                    axes1[i].set_xlabel('Peak Count', fontsize=12, fontweight='bold')
                    axes1[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements\nTop {len(traits)} Trait Associations', 
                                    fontsize=12, fontweight='bold', pad=20)
                    axes1[i].grid(True, alpha=0.3, linestyle='--', axis='x')
                    
                    # Set x-axis limits and ticks
                    max_count = max(counts)
                    axes1[i].set_xlim(0, max_count * 1.1)
                    axes1[i].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
                    
                    # Despine with offset
                    axes1[i].spines['top'].set_visible(False)
                    axes1[i].spines['right'].set_visible(False)
                    axes1[i].spines['left'].set_linewidth(1.5)
                    axes1[i].spines['bottom'].set_linewidth(1.5)
                    axes1[i].spines['left'].set_position(('outward', 10))
                    axes1[i].spines['bottom'].set_position(('outward', 10))
                    
                else:
                    axes1[i].text(0.5, 0.5, 'No trait data', ha='center', va='center', 
                               transform=axes1[i].transAxes, fontsize=14, fontweight='bold')
                    axes1[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements', 
                                    fontsize=12, fontweight='bold')
            else:
                axes1[i].text(0.5, 0.5, 'No intersections', ha='center', va='center', 
                           transform=axes1[i].transAxes, fontsize=14, fontweight='bold')
                axes1[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements', 
                                fontsize=12, fontweight='bold')
            
            # Despine with offset
            axes1[i].spines['top'].set_visible(False)
            axes1[i].spines['right'].set_visible(False)
            axes1[i].spines['left'].set_linewidth(1.5)
            axes1[i].spines['bottom'].set_linewidth(1.5)
            axes1[i].spines['left'].set_position(('outward', 10))
            axes1[i].spines['bottom'].set_position(('outward', 10))
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/trait_analysis_top10_combined.png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(f"{output_dir}/plots/trait_analysis_top10_combined.pdf", bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        # Figure 2: Top 35 traits (2x2 subplot)
        fig2, axes2 = plt.subplots(2, 2, figsize=(20, 12), facecolor='white')
        axes2 = axes2.ravel()
        
        for i, reg_type in enumerate(regulatory_files):
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                data = intersection_results[reg_type]
                trait_counts = analyze_traits(data)
                
                if trait_counts:
                    # Use prioritized traits (max 35 for combined plot)
                    prioritized_traits = prioritize_traits_for_plotting(trait_counts, max_traits=35)
                    
                    # Sort by count (descending order)
                    prioritized_traits = sorted(prioritized_traits, key=lambda x: x[1], reverse=True)
                    
                    traits, counts = zip(*prioritized_traits)
                    
                    # Clean trait names for better display
                    clean_traits = [clean_trait_name(trait) for trait in traits]
                    
                    # Ensure counts are integers
                    counts = [int(count) for count in counts]
                    
                    # Create horizontal bar plot with consistent bar width
                    y_positions = np.arange(len(clean_traits))
                    axes2[i].barh(y_positions, counts, height=0.6, alpha=0.8, 
                                 color=color_map[reg_type], edgecolor='black', linewidth=0.5)
                    axes2[i].set_yticks(y_positions)
                    
                    # Create y-axis labels with highlighting for non-ukb traits
                    y_labels = []
                    for j, (trait, clean_trait) in enumerate(zip(traits, clean_traits)):
                        if is_non_ukb_trait(trait):
                            y_labels.append(f"★ {clean_trait}")
                        else:
                            y_labels.append(clean_trait)
                    
                    axes2[i].set_yticklabels(y_labels, fontsize=8)  # Smaller font for more traits
                    
                    # Color non-ukb trait labels red
                    for j, (trait, label) in enumerate(zip(traits, y_labels)):
                        if is_non_ukb_trait(trait):
                            axes2[i].get_yticklabels()[j].set_color('red')
                            axes2[i].get_yticklabels()[j].set_fontweight('bold')
                    axes2[i].set_xlabel('Peak Count', fontsize=12, fontweight='bold')
                    axes2[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements\nTop {len(traits)} Trait Associations', 
                                    fontsize=12, fontweight='bold', pad=20)
                    axes2[i].grid(True, alpha=0.3, linestyle='--', axis='x')
                    
                    # Set x-axis limits and ticks
                    max_count = max(counts)
                    axes2[i].set_xlim(0, max_count * 1.1)
                    axes2[i].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
                    
                    # Despine with offset
                    axes2[i].spines['top'].set_visible(False)
                    axes2[i].spines['right'].set_visible(False)
                    axes2[i].spines['left'].set_linewidth(1.5)
                    axes2[i].spines['bottom'].set_linewidth(1.5)
                    axes2[i].spines['left'].set_position(('outward', 10))
                    axes2[i].spines['bottom'].set_position(('outward', 10))
                    
                else:
                    axes2[i].text(0.5, 0.5, 'No trait data', ha='center', va='center', 
                               transform=axes2[i].transAxes, fontsize=14, fontweight='bold')
                    axes2[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements', 
                                    fontsize=12, fontweight='bold')
            else:
                axes2[i].text(0.5, 0.5, 'No overlaps', ha='center', va='center', 
                           transform=axes2[i].transAxes, fontsize=14, fontweight='bold')
                axes2[i].set_title(f'{reg_type.replace("_", " ").title()} Regulatory Elements', 
                                fontsize=12, fontweight='bold')
            
            # Despine with offset
            axes2[i].spines['top'].set_visible(False)
            axes2[i].spines['right'].set_visible(False)
            axes2[i].spines['left'].set_linewidth(1.5)
            axes2[i].spines['bottom'].set_linewidth(1.5)
            axes2[i].spines['left'].set_position(('outward', 10))
            axes2[i].spines['bottom'].set_position(('outward', 10))
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/trait_analysis_top35_combined.png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(f"{output_dir}/plots/trait_analysis_top35_combined.pdf", bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
    except Exception as e:
        print(f"  Error creating combined trait plots: {e}")
        raise

def create_visualizations(intersection_stats, intersection_results, cwas_data, output_dir):
    """Create comprehensive visualizations"""
    print("\nSTEP 4: Creating visualizations...")
    
    try:
        # 1. Intersection Overview
        plot_intersection_overview(intersection_stats, cwas_data, output_dir)
        
        # 2. Pleiotropy Analysis
        plot_pleiotropy_analysis(intersection_results, output_dir)
        
        # 3. Trait Association Analysis
        plot_trait_analysis(intersection_results, output_dir)
        
        print("All visualizations saved to plots/ directory")
        
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        raise

def create_detailed_trait_analysis_file(intersection_results, cwas_data, output_dir):
    """Create detailed file with all trait and cQTL information for each regulatory region separately"""
    print("  Creating detailed trait analysis files for each regulatory region...")
    
    try:
        all_detailed_data = []
        
        # Create separate folder for detailed stats
        detailed_stats_dir = f"{output_dir}/detailed_stats"
        os.makedirs(detailed_stats_dir, exist_ok=True)
        print(f"  Created detailed stats directory: {detailed_stats_dir}")
        
        # Process all regulatory regions, even if they have no intersections
        all_regulatory_types = ['adult_specific', 'developmental_specific', 'fetal_specific', 'stem_cell_specific']
        
        # Debug: Show what regulatory types are available
        print(f"    Available regulatory types in intersection_results: {list(intersection_results.keys())}")
        for reg_type, data in intersection_results.items():
            print(f"      {reg_type}: {len(data)} intersections")
        
        for reg_type in all_regulatory_types:
            print(f"    Processing {reg_type} regulatory region...")
            
            # Check if this regulatory type has intersection data
            if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                results = intersection_results[reg_type]
                print(f"      Found {len(results)} intersections for {reg_type}")
                # Add regulatory type information
                results_copy = results.copy()
                results_copy['regulatory_type'] = reg_type
            else:
                # Create empty DataFrame for regulatory types with no intersections
                print(f"    No intersections found for {reg_type}, creating empty file...")
                results_copy = pd.DataFrame({'regulatory_type': [reg_type]})
                
            # Add trait analysis for each row
            trait_analysis = []
            if len(results_copy) > 0:
                # Check if all_descriptions column exists, if not create it
                if 'all_descriptions' not in results_copy.columns:
                    results_copy['all_descriptions'] = ''
                    print(f"    Warning: 'all_descriptions' column not found for {reg_type}, creating empty column")
                
                for idx, row in results_copy.iterrows():
                    if pd.notna(row.get('all_descriptions')) and row.get('all_descriptions') != '':
                        traits = row['all_descriptions'].split('; ')
                        
                        trait_analysis.append({
                            'total_traits': len(traits),
                            'all_traits': '; '.join(traits)
                        })
                    else:
                        trait_analysis.append({
                            'total_traits': 0,
                            'all_traits': ''
                        })
            else:
                # For empty DataFrames, add empty trait analysis
                trait_analysis.append({
                    'total_traits': 0,
                    'all_traits': ''
                })
            
            # Add required columns for empty DataFrames
            if len(results_copy) == 0:
                results_copy['pleiotropy'] = ['trait-specific']  # Default value
                results_copy['all_descriptions'] = ['']
            
            # Add trait analysis columns
            trait_df = pd.DataFrame(trait_analysis)
            results_with_traits = pd.concat([results_copy, trait_df], axis=1)
            
            # Save individual regulatory region detailed file (CSV)
            individual_file = f"{detailed_stats_dir}/detailed_trait_analysis_{reg_type}.csv"
            results_with_traits.to_csv(individual_file, index=False)
            
            # Save individual regulatory region detailed file (TSV)
            individual_tsv_file = f"{detailed_stats_dir}/detailed_trait_analysis_{reg_type}.tsv"
            results_with_traits.to_csv(individual_tsv_file, sep='\t', index=False)
            
            print(f"    Saved detailed analysis for {reg_type}: {individual_file}")
            print(f"    Saved detailed TSV for {reg_type}: {individual_tsv_file}")
            print(f"    Records for {reg_type}: {len(results_with_traits)}")
            
            # Create trait statistics for this regulatory region
            trait_counts = analyze_traits(results_with_traits)
            
            # Get plotted traits (prioritized and sorted)
            prioritized_traits = prioritize_traits_for_plotting(trait_counts, max_traits=35)
            prioritized_traits = sorted(prioritized_traits, key=lambda x: x[1], reverse=True)
            
            # Create plotted traits information
            plotted_traits_info = []
            for trait, count in prioritized_traits:
                clean_trait = clean_trait_name(trait)
                plotted_traits_info.append({
                    'original_trait': str(trait),
                    'clean_trait': str(clean_trait),
                    'count': int(count)
                })
            
            # Convert numpy types to Python native types for JSON serialization
            # Handle pleiotropy distribution safely
            pleiotropy_dist = {}
            if 'pleiotropy' in results_with_traits.columns:
                pleiotropy_dist = {str(k): int(v) for k, v in results_with_traits['pleiotropy'].value_counts().to_dict().items()}
            else:
                pleiotropy_dist = {'trait-specific': 0, 'pleiotropic': 0}
            
            trait_stats = {
                'regulatory_type': reg_type,
                'total_intersections': int(len(results_with_traits)),
                'total_unique_traits': int(len(trait_counts)),
                'plotted_traits_count': int(len(plotted_traits_info)),
                'all_traits_count': int(len(trait_counts)),
                'trait_frequency': {str(k): int(v) for k, v in dict(trait_counts).items()},
                'plotted_traits': plotted_traits_info,
                'pleiotropy_distribution': pleiotropy_dist,
                'total_traits_analyzed': int(results_with_traits['total_traits'].sum()) if 'total_traits' in results_with_traits.columns else 0
            }
            
            # Save individual regulatory region statistics
            stats_file = f"{detailed_stats_dir}/trait_analysis_stats_{reg_type}.json"
            import json
            with open(stats_file, 'w') as f:
                json.dump(trait_stats, f, indent=2)
            
            # Save plotted traits as CSV and TSV for easier access (top traits only)
            plotted_traits_df = pd.DataFrame(plotted_traits_info)
            plotted_traits_file = f"{detailed_stats_dir}/plotted_traits_{reg_type}.csv"
            plotted_traits_df.to_csv(plotted_traits_file, index=False)
            
            plotted_traits_tsv_file = f"{detailed_stats_dir}/plotted_traits_{reg_type}.tsv"
            plotted_traits_df.to_csv(plotted_traits_tsv_file, sep='\t', index=False)
            
            # Create ALL traits file (not just the plotted ones)
            all_traits_info = []
            for trait, count in trait_counts:
                clean_trait = clean_trait_name(trait)
                all_traits_info.append({
                    'original_trait': str(trait),
                    'clean_trait': str(clean_trait),
                    'count': int(count),
                    'included_in_plot': trait in [t[0] for t in prioritized_traits]
                })
            
            # Sort by count (descending order)
            all_traits_info = sorted(all_traits_info, key=lambda x: x['count'], reverse=True)
            
            # Save ALL traits as CSV and TSV
            all_traits_df = pd.DataFrame(all_traits_info)
            all_traits_file = f"{detailed_stats_dir}/all_traits_{reg_type}.csv"
            all_traits_df.to_csv(all_traits_file, index=False)
            
            all_traits_tsv_file = f"{detailed_stats_dir}/all_traits_{reg_type}.tsv"
            all_traits_df.to_csv(all_traits_tsv_file, sep='\t', index=False)
            
            print(f"    Saved statistics for {reg_type}: {stats_file}")
            print(f"    Saved plotted traits for {reg_type}: {plotted_traits_file}")
            print(f"    Saved plotted traits TSV for {reg_type}: {plotted_traits_tsv_file}")
            print(f"    Saved ALL traits for {reg_type}: {all_traits_file}")
            print(f"    Saved ALL traits TSV for {reg_type}: {all_traits_tsv_file}")
            
            all_detailed_data.append(results_with_traits)
        
        if all_detailed_data:
            # Combine all regulatory types for overall analysis
            combined_detailed = pd.concat(all_detailed_data, ignore_index=True)
            
            # Save combined detailed file (CSV)
            combined_file = f"{output_dir}/detailed_trait_analysis_combined.csv"
            combined_detailed.to_csv(combined_file, index=False)
            
            # Save combined detailed file (TSV)
            combined_tsv_file = f"{output_dir}/detailed_trait_analysis_combined.tsv"
            combined_detailed.to_csv(combined_tsv_file, sep='\t', index=False)
            
            print(f"  Combined detailed trait analysis saved to: {combined_file}")
            print(f"  Combined detailed TSV saved to: {combined_tsv_file}")
            print(f"  Total detailed records: {len(combined_detailed)}")
            
            # Create overall summary statistics (simplified - no cancer/non-ukb separation)
            overall_summary = {
                'total_intersections': int(len(combined_detailed)),
                'regulatory_types': {str(k): int(v) for k, v in combined_detailed['regulatory_type'].value_counts().to_dict().items()},
                'pleiotropy_distribution': {str(k): int(v) for k, v in combined_detailed['pleiotropy'].value_counts().to_dict().items()},
                'total_traits_analyzed': int(combined_detailed['total_traits'].sum()),
                'total_unique_traits': len(combined_detailed['all_traits'].str.split('; ').explode().unique()) if len(combined_detailed) > 0 else 0
            }
            
            # Save overall summary statistics
            overall_summary_file = f"{output_dir}/trait_analysis_summary_overall.json"
            with open(overall_summary_file, 'w') as f:
                json.dump(overall_summary, f, indent=2)
            
            print(f"  Overall trait analysis summary saved to: {overall_summary_file}")
            
            return combined_detailed, overall_summary
        else:
            print("  No intersection data available for detailed analysis")
            return pd.DataFrame(), {}
            
    except Exception as e:
        print(f"  Error creating detailed trait analysis file: {e}")
        raise

def generate_summary_report(intersection_stats, intersection_results, cwas_data, output_dir):
    """Generate a comprehensive summary report"""
    print("\nSTEP 5: Generating summary report...")
    
    try:
        report_file = f"{output_dir}/cfCWAS_Regulatory_Analysis_Report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("cfCWAS REGULATORY ELEMENT INTERSECTION ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("ANALYSIS OVERVIEW\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total cfCWAS loci analyzed: {len(cwas_data)}\n")
            f.write("Specific regulatory element types: adult_specific, developmental_specific, fetal_specific, stem_cell_specific\n")
            f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("Tool used: pybedtools\n\n")
            
            f.write("PLEIOTROPY SUMMARY\n")
            f.write("-" * 40 + "\n")
            pleiotropy_summary = cwas_data['pleiotropy'].value_counts()
            for pleiotropy_type, count in pleiotropy_summary.items():
                percentage = count / len(cwas_data) * 100
                f.write(f"{pleiotropy_type}: {count} ({percentage:.1f}%)\n")
            f.write("\n")
            
            f.write("INTERSECTION RESULTS\n")
            f.write("-" * 40 + "\n")
            for reg_type, stats in intersection_stats.items():
                f.write(f"\n{reg_type.replace('_', ' ').upper()} REGULATORY ELEMENTS:\n")
                f.write(f"  Total intersections: {stats['total_intersections']}\n")
                f.write(f"  Unique cfCWAS loci: {stats['unique_cwas_loci']}\n")
                f.write(f"  Unique cQTL elements: {stats['unique_cQTL_elements']}\n")
                if stats['unique_cwas_loci'] > 0:
                    intersection_rate = stats['unique_cwas_loci'] / len(cwas_data) * 100
                    f.write(f"  Intersection rate: {intersection_rate:.1f}%\n")
                
                # Pleiotropy analysis for this regulatory type
                if reg_type in intersection_results and len(intersection_results[reg_type]) > 0:
                    intersected_data = intersection_results[reg_type]
                    pleiotropy_counts = intersected_data['pleiotropy'].value_counts()
                    f.write("  Pleiotropy distribution:\n")
                    for pleiotropy_type, count in pleiotropy_counts.items():
                        percentage = count / stats['unique_cwas_loci'] * 100
                        f.write(f"    {pleiotropy_type}: {count} ({percentage:.1f}%)\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("END OF REPORT\n")
            f.write("=" * 80 + "\n")
        
        print(f"Summary report saved to: {report_file}")
        
    except Exception as e:
        print(f"Error generating summary report: {e}")
        raise

def main():
    """Main function to run all analysis configurations serially"""
    
    # Define all analysis configurations
    # Update these paths to your data directories
    # Example configurations (uncomment and modify as needed):
    #analysis_configs = [
    #    {
    #        'name': 'cQTL_cfCWAS_UKB',
    #        'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/cfcQTLs_unique_to_WBCcQTLs",
    #        'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.UKB.pleiotropy.csv",
    #        'output_dir': "/path/to/output/NewTest_cQTL_cfCWAS_H3K4ME3_regulatory_intersect_UKBB"
    #    },
    #    {
    #        'name': 'cQTL_cfCWAS_nonUKB',
    #        'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/cfcQTLs_unique_to_WBCcQTLs",
    #        'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.non-UKB.pleiotropy.csv",
    #        'output_dir': "/path/to/output/NewTest_cQTL_cfCWAS_H3K4ME3_regulatory_intersect"
    #    },
    #    {
    #        'name': 'WBCcQTL_cfCWAS_nonUKB',
    #        'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/WBCcQTLs_unique_to_cfcQTLs",
    #        'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.non-UKB.pleiotropy.csv",
    #        'output_dir': "/path/to/output/NewTest_WBCcQTL_cfCWAS_H3K4ME3_regulatory_intersect"
    #    },
    #    {
    #        'name': 'WBCcQTL_cfCWAS_UKB',
    #        'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/WBCcQTLs_unique_to_cfcQTLs",
    #        'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.UKB.pleiotropy.csv",
    #        'output_dir': "/path/to/output/NewTest_WBCcQTL_cfCWAS_H3K4ME3_regulatory_intersect_UKBB"
    #    }
    #]
    
    # Example configuration for combined analysis
    # Update paths to your actual data directories
    analysis_configs = [
        {
            'name': 'cQTL_cfCWAS_combined',
            'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/cfcQTLs_unique_to_WBCcQTLs",
            'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.all.pleiotropy.csv",
            'output_dir': "/path/to/output/NewTest_cQTL_cfCWAS_H3K4ME3_regulatory_intersect_COMBINED"
        },
        {
            'name': 'WBCcQTL_cfCWAS_combined',
            'regulatory_dir': "/path/to/workspace/NewTest_cQTL_H3K4me3_qvalBased_V2/WBCcQTLs_unique_to_cfcQTLs",
            'cwas_file': "/path/to/data/CWAS/cfChIP_H3K4me3.cwas.peaks.all.pleiotropy.csv",
            'output_dir': "/path/to/output/NewTest_WBCcQTL_cfCWAS_H3K4ME3_regulatory_intersect_COMBINED"
        }
    ]
    
    print("cfCWAS Regulatory Element Intersection Analysis - Automated Serial Execution")
    print("=" * 80)
    print(f"Total configurations to run: {len(analysis_configs)}")
    print("=" * 80)
    
    # Store results for all configurations
    all_results = {}
    
    for i, config in enumerate(analysis_configs, 1):
        print(f"\n{'='*20} ANALYSIS {i}/{len(analysis_configs)}: {config['name']} {'='*20}")
        print(f"Regulatory Directory: {config['regulatory_dir']}")
        print(f"cfCWAS File: {config['cwas_file']}")
        print(f"Output Directory: {config['output_dir']}")
        print("=" * 80)
        
        # Check if files exist
        if not os.path.exists(config['cwas_file']):
            print(f"Error: cfCWAS file not found: {config['cwas_file']}")
            print(f"Skipping analysis {i}: {config['name']}")
            all_results[config['name']] = {'status': 'failed', 'error': 'cfCWAS file not found'}
            continue
        
        if not os.path.exists(config['regulatory_dir']):
            print(f"Error: Regulatory directory not found: {config['regulatory_dir']}")
            print(f"Skipping analysis {i}: {config['name']}")
            all_results[config['name']] = {'status': 'failed', 'error': 'Regulatory directory not found'}
            continue
        
        try:
            # Create output directory structure
            os.makedirs(config['output_dir'], exist_ok=True)
            os.makedirs(f"{config['output_dir']}/intersections", exist_ok=True)
            os.makedirs(f"{config['output_dir']}/plots", exist_ok=True)
            os.makedirs(f"{config['output_dir']}/statistics", exist_ok=True)
            print(f"Created output directory: {config['output_dir']}")
            
            # STEP 1: Load cfCWAS data
            cwas_data, cwas_bed_file = load_cwas_data(config['cwas_file'], config['output_dir'])
            
            # STEP 2: Perform intersections
            intersection_results, intersection_stats = perform_intersections(
                cwas_bed_file, config['regulatory_dir'], config['output_dir'], cwas_data
            )
            
            # STEP 3: Generate statistics
            stats_df = generate_statistics(intersection_stats, cwas_data, intersection_results, config['output_dir'])
            
            # STEP 4: Create visualizations
            create_visualizations(intersection_stats, intersection_results, cwas_data, config['output_dir'])
            
            # STEP 5: Create detailed trait analysis file
            detailed_trait_data, trait_summary = create_detailed_trait_analysis_file(
                intersection_results, cwas_data, config['output_dir']
            )
            
            # STEP 6: Generate summary report
            generate_summary_report(intersection_stats, intersection_results, cwas_data, config['output_dir'])
            
            print(f"\n{'='*20} ANALYSIS {i} COMPLETE: {config['name']} {'='*20}")
            print(f"Results saved to: {config['output_dir']}")
            print(f"Total CWAS loci: {len(cwas_data)}")
            print(f"Regulatory types analyzed: {len(intersection_stats)}")
            
            # Store results
            all_results[config['name']] = {
                'status': 'success',
                'stats_df': stats_df,
                'intersection_results': intersection_results,
                'intersection_stats': intersection_stats,
                'cwas_data': cwas_data,
                'output_dir': config['output_dir']
            }
            
        except Exception as e:
            print(f"\nANALYSIS {i} FAILED: {config['name']}")
            print(f"Error: {e}")
            print("Please check the error messages above and fix any issues.")
            all_results[config['name']] = {'status': 'failed', 'error': str(e)}
            continue
    
    # Print final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY - ALL ANALYSES")
    print("=" * 80)
    
    successful_analyses = 0
    failed_analyses = 0
    
    for config_name, result in all_results.items():
        if result['status'] == 'success':
            successful_analyses += 1
            print(f"{config_name}: SUCCESS")
        else:
            failed_analyses += 1
            print(f"{config_name}: FAILED - {result.get('error', 'Unknown error')}")
    
    print(f"\nTotal Analyses: {len(analysis_configs)}")
    print(f"Successful: {successful_analyses}")
    print(f"Failed: {failed_analyses}")
    print("=" * 80)
    
    return all_results

if __name__ == "__main__":
    all_results = main() 


