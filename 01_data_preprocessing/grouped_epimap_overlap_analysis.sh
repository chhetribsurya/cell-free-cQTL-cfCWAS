#!/bin/bash

module load gcc/9.2.0
module load bedtools/2.30.0

# Update this path to your working directory
# cd /path/to/your/workspace
mkdir -p work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/
# Extract unique groups from the third column of the CSV file
cut -d ',' -f3 data/EpiMap/EpiMap_info.csv | sort | uniq | while read -r line; do
    # Sanitize the group name
    sanitized=$(echo "$line" | tr -cd '[:alnum:]_')
    # Extract first column for rows matching the group and save to a file
    awk -F',' -v val="$line" '$3 == val' data/EpiMap/EpiMap_info.csv | cut -d ',' -f1 > "${sanitized}.txt"
done

# Process each group
ls *.txt | sed 's/.txt//g' | while read -r group; do
    # Concatenate enhancer files for each group into a temporary file
    while read -r line; do
        cat "data/EpiMap/epigenome_18/EnhA1_EnhA2_EnhG1_EnhG2/1.0/H3K27ac_${line}_18_core_K27ac_dense.bed"
    done < "${group}.txt" > "${group}.tmp"

    # Merge peaks with a threshold
    threshold=1
    sh scripts/merge.peaks.sh "${group}.tmp" $threshold "${group}.merged.bed"
done

mv *.merged.bed  work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/
rm *.tmp
rm *.txt
rm counts.bed

> work/EpiMap_overlap/EnhA1_EnhA2_EnhG1_EnhG2_grouped_output_ratios.txt
# Calculate ratios for cfChIP_H3K27ac data
ls work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/ | sed 's/.merged.bed//g' | while read -r line; do
    count=$(sed 's/^/chr/g' data/cfChIP_H3K27ac/combined.sig.bed | bedtools intersect -a - -b work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/${line}.merged.bed -wa | sort -k1,1n -k2,2n | uniq | wc -l)
    sum=$(wc -l < data/cfChIP_H3K27ac/combined.sig.bed)
    ratio=$(echo "scale=4; $count/$sum" | bc)
    echo -e "$line\t$ratio\tcfChIP_H3K27ac" >> work/EpiMap_overlap/EnhA1_EnhA2_EnhG1_EnhG2_grouped_output_ratios.txt
done

# Calculate ratios for whole-blood data
ls work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/ | sed 's/.merged.bed//g' | while read -r line; do
    count=$(sed 's/^/chr/g' data/whole-blood/combined.sig.bed | bedtools intersect -a - -b work/EpiMap_overlap/grouped_merged_EpiMap_EnhA1_EnhA2_EnhG1_EnhG2/${line}.merged.bed -wa | sort -k1,1n -k2,2n | uniq | wc -l)
    sum=$(wc -l < data/whole-blood/combined.sig.bed)
    ratio=$(echo "scale=4; $count/$sum" | bc)
    echo -e "$line\t$ratio\twhole-blood" >> work/EpiMap_overlap/EnhA1_EnhA2_EnhG1_EnhG2_grouped_output_ratios.txt
done

