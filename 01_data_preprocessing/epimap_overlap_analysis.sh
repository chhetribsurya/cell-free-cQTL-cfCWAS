#!/bin/bash

# Update this path to your working directory
# cd /path/to/your/workspace

# Load necessary modules
module load gcc/9.2.0
module load bedtools/2.30.0

# Define directories and files
# Input arguments
raw_dir="data/EpiMap/epigenome_18"        # First argument: raw directory
patterns=$1       # Second argument: patterns (comma-separated)

# Create output directory name based on patterns
out=$(echo "$patterns" | tr ',' '_')  # Replace commas with underscores
output_dir_states="${out}/1.0"

# Create necessary directories
mkdir -p "${raw_dir}/${output_dir_states}"

# Define other directories and files
enhancer_dir="${raw_dir}/${output_dir_states}"
# Update these paths to your data files
consensus_peaks_cf="/path/to/your/data/cfChIP_H3K27ac/analysis/consensus_peaks/consensus.peaks.stratas.tsv"
consensus_peaks_wb="/path/to/your/data/whole-blood/cwas/analysis/consensus_peaks/consensus.peaks.stratas.tsv"
label_file="data/EpiMap/label.csv"
output_dir="work/EpiMap_overlap"
tmp_dir="work/tmp"
combined_output="work/EpiMap_overlap/${out}_combined_results.tsv"

# Log directories and settings
echo "Raw directory: $raw_dir"
echo "Patterns: $patterns"
echo "Output directory states: $output_dir_states"
echo "Enhancer directory: $enhancer_dir"
echo "Combined output file: $combined_output"

# Create necessary directories
mkdir -p "$output_dir/cf_k27"
mkdir -p "$output_dir/wb_k27"
mkdir -p "$tmp_dir"

# Temporary output files
line_tmp="$tmp_dir/line.tmp"
name_tmp="$tmp_dir/name.tmp"
wb_count_tmp="$tmp_dir/wb.count.tmp"
cf_count_tmp="$tmp_dir/cf.count.tmp"
total_count_tmp="$tmp_dir/total.count.tmp"
downsampled_output="work/EpiMap_overlap/whole_blood_downsampled_by_chrom.bed"

# Clear temporary files
> "$line_tmp"
> "$name_tmp"
> "$wb_count_tmp"
> "$cf_count_tmp"
> "$total_count_tmp"
> "$combined_output"

# Step 1: Select chromatin states and merge based on patterns

# Convert patterns from a comma-separated string into an array
IFS=',' read -r -a pattern_array <<< "$patterns"

for file in ${raw_dir}/*.gz; do
    # Extract the base name of the file
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed.gz//g')
    output_file="${raw_dir}/${output_dir_states}/H3K27ac_${base_name}_18_core_K27ac_dense.bed"

    # Process the file for each pattern in the array
    for pattern in "${pattern_array[@]}"; do
        zcat "$file" | grep "$pattern" | sort -k1,1 -k2,2n | bedtools merge -i - >> "$output_file"
    done

    echo "Processed: $file -> $output_file"
done


# Step 2: Find overlap between cf peaks and chromatin states
for file in "$enhancer_dir"/*_18_core_K27ac_dense.bed; do
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')
    filtered_consensus="${output_dir}/cf_k27/filtered_consensus_peaks.bed"
    sed 's/^/chr/g' "$consensus_peaks_cf" | cut -f1-3 | grep -v "CHR" > "$filtered_consensus"
    
    output_file="${output_dir}/cf_k27/cfChIP_${base_name}_${out}_1.0_overlap.bed"
    bedtools intersect -a "$file" -b  "$filtered_consensus" | uniq > "$output_file"
    echo "Processed: $file -> $output_file"
done

rm "${output_dir}/cf_k27/filtered_consensus_peaks.bed"

# Step 3: Downsample whole-blood peaks by chromosome counts
chromosomes=$(cut -f1 "$consensus_peaks_cf" | grep -v -e "X" -e "Y" -e "M" -e "CHR" | sort -k1,1n | uniq)
cut -f1 "$consensus_peaks_cf" | grep -v -e "X" -e "Y" -e "M" -e "CHR" | sort -k1,1n | uniq -c | awk '{print $2"\t"$1}' > "$tmp_dir/target_chromosome_counts.txt"

> "$downsampled_output"
while read -r chrom count; do
    grep -w "^$chrom" "$consensus_peaks_wb" > "$tmp_dir/temp_chrom.bed"
    shuf "$tmp_dir/temp_chrom.bed" | head -n "$count" >> "$downsampled_output"
done < "$tmp_dir/target_chromosome_counts.txt"
rm "$tmp_dir/temp_chrom.bed"

# Find overlaps for whole-blood peaks
for file in "$enhancer_dir"/*_18_core_K27ac_dense.bed; do
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')
    filtered_consensus="${output_dir}/wb_k27/filtered_consensus_peaks.bed"
    sed 's/^/chr/g' "$downsampled_output" | cut -f1-3 | grep -v "CHR" > "$filtered_consensus"
    
    output_file="${output_dir}/wb_k27/wbc_${base_name}_${out}_1.0_overlap.bed"
    bedtools intersect -a "$file" -b  "$filtered_consensus" |  uniq > "$output_file"
    echo "Processed: $file -> $output_file"
done

# Clean up intermediate files
rm "$downsampled_output"

# Step 4: Combine results into a single file
ls "$enhancer_dir" | sed 's/_18_core_K27ac_dense.bed//g' | sed 's/H3K27ac_//g' | while read -r line; do
    echo "$line" >> "$line_tmp"
    grep -w "$line" "$label_file" | cut -f2 >> "$name_tmp"
    wc -l < "$output_dir/wb_k27/wbc_H3K27ac_${line}_${out}_1.0_overlap.bed" >> "$wb_count_tmp"
    wc -l < "$output_dir/cf_k27/cfChIP_H3K27ac_${line}_${out}_1.0_overlap.bed" >> "$cf_count_tmp"
    wc -l < "$enhancer_dir/H3K27ac_${line}_18_core_K27ac_dense.bed" >> "$total_count_tmp"
done

# Combine all results with a new column for total counts
paste "$line_tmp" "$wb_count_tmp" "$cf_count_tmp" "$total_count_tmp" "$name_tmp" | \
    awk 'BEGIN {print "EID\tWB_count\tCF_count\tTotal_count\tNAME"} {print $0}' > "$combined_output"

echo "Combined results saved to: $combined_output"

