#!/usr/bin/env bash
module load gcc/9.2.0
module load bedtools/2.30.0

#generate all the chromatin states and downsampled peak:
# Update this path to your working directory
# cd /path/to/your/workspace

query=$1 #e.g., data/cfChIP_H3K27ac/consensus.peaks.bed
assay=$2 #cfChIP.H3K27ac

# 1) Ensure query is sorted
sort -k1,1 -k2,2n $query | cut -f 1-3> query.sorted.bed

# 2) Make a place for per‐ref cols
mkdir -p tmp_counts

# 3) Loop over your 98 refs (adjust the pattern to point at them)
for f in data/EpiMap/epigenome_18/EnhA1_EnhA2_TssA_EnhG1_EnhG2/1.0/*.noChr.sorted.merged.bed; do
  # strip off the directory + extension
  name=$(basename "$f" .noChr.sorted.merged.bed)
  # get a 0/1 vector of whether each query peak overlaps this ref
  bedtools intersect -c \
    -a query.sorted.bed \
    -b "$f" \
  | cut -f4 \
  | sed 's/^[1-9]/1/' \
  > tmp_counts/"${name}".col
done

# 4) Paste coords + all 98 cols, then sum across them:
(printf "chr\tstart\tend\t" && basename -a tmp_counts/*.col   | sed 's/^H3K27ac_//'   | sed 's/_18_core_K27ac_dense.col$//'   | paste -sd '\t') > header.txt
paste <(cut -f1-3 query.sorted.bed) tmp_counts/*.col > body.txt
cat header.txt body.txt > merged_output.tsv

paste \
  <(cut -f1-3 query.sorted.bed) \
  tmp_counts/*.col \
| awk '{
    sum=0
    for(i=4;i<=NF;i++) sum+=($i>0)
    print $1"\t"$2"\t"$3"\t"sum
  }' \
> query.detect_counts.bed

mv merged_output.tsv work/EpiMap_overlap/${assay}.peak.overlaped.tsv
mv query.detect_counts.bed work/EpiMap_overlap/${assay}.peak.freq.bed

rm header.txt body.txt
rm -r tmp_counts

echo 'Done: results in .peak.overlaped.tsv and .peak.freq.bed'