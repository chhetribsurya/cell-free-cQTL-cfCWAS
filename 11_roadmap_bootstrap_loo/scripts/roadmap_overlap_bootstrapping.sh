#!/bin/bash
set -euo pipefail

cd /n/scratch/users/z/ziz597/cwas/summary/work/revision/epimap_overlap


#module load gcc/9.2.0
module load bedtools/2.31.0

# -------------------------
# Inputs
# -------------------------
raw_dir="/n/scratch/users/z/ziz597/cwas/summary/data/EpiMap/epigenome_18"
patterns=$1                   # comma-separated (e.g. "EnhA1,EnhA2,EnhG1,EnhG2")
B=${2:-200}                   # bootstrap replicates
SEED0=${3:-1}                 # base seed
FRAC=${4:-0.8}  



out=$(echo "$patterns" | tr ',' '_')
output_dir_states="${out}/1.0"

enhancer_dir="${raw_dir}/${output_dir_states}"

consensus_peaks_cf="../enrichment/data/background/cfChIP_H3K27ac.consensus.peaks.bed"
consensus_peaks_wb="../enrichment/data/background/wb.consensus.peaks.bed"

label_file="/n/scratch/users/z/ziz597/cwas/summary/data/EpiMap/label.csv"

output_dir="work/"
tmp_dir="work/tmp/${out}"
mkdir -p "${raw_dir}/${output_dir_states}" "$output_dir" "$tmp_dir"

combined_output="${output_dir}/${out}_combined_results.bootstrap.tsv"
boot_long="${tmp_dir}/${out}.boot.long.tsv"

echo "============================================================"
echo "patterns: $patterns"
echo "B: $B"
echo "SEED0: $SEED0"
echo "enhancer_dir: $enhancer_dir"
echo "combined_output: $combined_output"
echo "============================================================"

# -------------------------
# Helpers: reproducible shuf random-source
# -------------------------
rand_src() {
  local seed="$1"
  # Deterministic stream of bytes from seed (works well on clusters)
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

# -------------------------
# Step 1: Build enhancer BEDs if missing
# -------------------------
IFS=',' read -r -a pattern_array <<< "$patterns"

if ls "${enhancer_dir}"/*_18_core_K27ac_dense.bed >/dev/null 2>&1; then
  echo "[Step1] Enhancer BEDs exist, skipping build."
else
  echo "[Step1] Building enhancer BEDs..."
  mkdir -p "${enhancer_dir}"
  for file in ${raw_dir}/*.gz; do
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed.gz//g')
    output_file="${raw_dir}/${output_dir_states}/H3K27ac_${base_name}_18_core_K27ac_dense.bed"
    > "$output_file"

    for pattern in "${pattern_array[@]}"; do
      zcat "$file" | grep "$pattern" | sort -k1,1 -k2,2n | bedtools merge -i - >> "$output_file"
    done

    echo "Processed: $file -> $output_file"
  done
fi

enh_files=($(ls "${enhancer_dir}"/*_18_core_K27ac_dense.bed))
if [ ${#enh_files[@]} -eq 0 ]; then
  echo "[ERROR] No enhancer files found in ${enhancer_dir}"
  exit 1
fi

# -------------------------
# Step 2: Prepare per-chrom files for CF and WB + target chrom counts
# -------------------------
echo "[Step2] Pre-splitting CF and WB by chromosome..."

cf_dir="${tmp_dir}/cf_by_chrom"
wb_dir="${tmp_dir}/wb_by_chrom"
mkdir -p "$cf_dir" "$wb_dir"

# Remove headers + sex/MT; keep 3 columns for CF; keep full line for WB (but we only need first 3 later)
awk 'BEGIN{OFS="\t"} $1!="CHR" && $1!="X" && $1!="Y" && $1!="M"{print $1,$2,$3}' "$consensus_peaks_cf" \
  | awk -v D="$cf_dir" 'BEGIN{OFS="\t"}{f=D"/"$1".bed"; print $0 >> f}' 

awk 'BEGIN{OFS="\t"} $1!="CHR" && $1!="X" && $1!="Y" && $1!="M"{print $0}' "$consensus_peaks_wb" \
  | awk -v D="$wb_dir" 'BEGIN{OFS="\t"}{f=D"/"$1".tsv"; print $0 >> f}'

# target counts from CF
target_counts="${tmp_dir}/target_chromosome_counts.txt"

cut -f1 "$consensus_peaks_cf" \
  | grep -v -e "CHR" -e "X" -e "Y" -e "M" \
  | sort -k1,1n \
  | uniq -c \
  | awk -v frac="$FRAC" '
      BEGIN{OFS="\t"}
      {
        chrom=$2; n=$1;
        k=int(n*frac);
        if(k<1) k=1;
        print chrom,k
      }' \
  > "$target_counts"
  
chrom_list="${tmp_dir}/chrom_list.txt"
cut -f1 "$target_counts" > "$chrom_list"

echo "[Step2] Target chrom counts:"
head "$target_counts"

# -------------------------
# Step 3: Bootstrap loop
# Output long table: rep  EID  WB_count  CF_count
# -------------------------
echo "[Step3] Bootstrapping overlaps..."
> "$boot_long"

for b in $(seq 1 "$B"); do
  seed=$((SEED0 + b))

  cf_boot="${tmp_dir}/cf.boot.rep${b}.bed"       # nochr, 3-col
  wb_down="${tmp_dir}/wb.down.rep${b}.bed"       # nochr, 3-col
  > "$cf_boot"
  > "$wb_down"

  # --- Bootstrap CF within chrom (WITH replacement) ---
  while read -r chrom cnt; do
    cf_chrom_file="${cf_dir}/${chrom}.bed"
    if [ ! -s "$cf_chrom_file" ]; then
      continue
    fi

    # Sample with replacement using awk
    shuf --random-source=<(rand_src "${seed}_${chrom}_cf") "$cf_chrom_file" \
    | head -n "$cnt" \
    >> "$cf_boot"


  done < "$target_counts"

  # --- Downsample WB within chrom (WITHOUT replacement; seeded shuf) ---
  while read -r chrom cnt; do
    wb_chrom_file="${wb_dir}/${chrom}.tsv"
    if [ ! -s "$wb_chrom_file" ]; then
      continue
    fi

    # keep first 3 cols only
    shuf --random-source=<(rand_src "${seed}_${chrom}_wb") "$wb_chrom_file" \
      | head -n "$cnt" \
      | cut -f1-3 \
      >> "$wb_down"

  done < "$target_counts"

  # add chr prefix for bedtools intersect (enhancer beds are chr-prefixed)
  cf_boot_chr="${tmp_dir}/cf.boot.rep${b}.chr.bed"
  wb_down_chr="${tmp_dir}/wb.down.rep${b}.chr.bed"
  sed 's/^/chr/' "$cf_boot" > "$cf_boot_chr"
  sed 's/^/chr/' "$wb_down" > "$wb_down_chr"

  # --- Count overlaps for each enhancer file ---
  for file in "${enh_files[@]}"; do
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')

    wb_count=$(bedtools intersect -u -a "$file" -b "$wb_down_chr" | wc -l | awk '{print $1}')
    cf_count=$(bedtools intersect -u -a "$file" -b "$cf_boot_chr" | wc -l | awk '{print $1}')

    echo -e "${b}\t${base_name}\t${wb_count}\t${cf_count}" >> "$boot_long"
  done

  rm -f "$cf_boot" "$wb_down" "$cf_boot_chr" "$wb_down_chr"

  if (( b % 25 == 0 )); then
    echo "  finished bootstrap replicate $b / $B"
  fi
done

# -------------------------
# Step 4: deterministic Total_count + NAME
# -------------------------
echo "[Step4] Total_count + NAME..."

total_counts_tsv="${tmp_dir}/${out}.total_counts.tsv"
name_tsv="${tmp_dir}/${out}.name.tsv"
> "$total_counts_tsv"
> "$name_tsv"

for file in "${enh_files[@]}"; do
  base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')

  total_count=$(wc -l < "$file" | awk '{print $1}')
  echo -e "${base_name}\t${total_count}" >> "$total_counts_tsv"

  # same behavior as your original (label.csv appears CSV)
  # base_name like H3K27ac_<EID>
  eid=$(echo "$base_name" | sed 's/^H3K27ac_//')
  name=$(grep -w "$eid" "$label_file" | head -n 1 | cut -d',' -f2)
  echo -e "${base_name}\t${name}" >> "$name_tsv"
done

# -------------------------
# Step 5: summarize in R (no python)
# -------------------------
echo "[Step5] Summarizing bootstrap distribution in R..."

Rscript scripts/summarize_epimap_bootstrap.R \
  "$boot_long" \
  "$total_counts_tsv" \
  "$name_tsv" \
  "$combined_output"

echo "============================================================"
echo "DONE. Results: $combined_output"
echo "============================================================"

