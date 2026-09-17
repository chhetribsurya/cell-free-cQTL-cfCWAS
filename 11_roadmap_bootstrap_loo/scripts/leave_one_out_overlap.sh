#!/bin/bash
#===============================================================================
# Leave-one-cancer-type/site-out overlap analysis
#
# For each leave-out group (cancer type or site):
#   1. Build consensus peaks from remaining samples
#   2. Downsample WBC peaks to match cfChIP peak count per chromosome
#   3. Compute overlap with EpiMap enhancer states
#
# Usage: bash leave_one_out_overlap.sh <patterns> [seed]
#   patterns: comma-separated ChromHMM states (e.g. "EnhA1,EnhA2,EnhG1,EnhG2")
#   seed: random seed for WBC downsampling (default: 42)
#
# Requires: bedtools, Rscript, openpyxl (via python for xlsx reading)
#===============================================================================

set -euo pipefail

cd /n/scratch/users/z/ziz597/cwas/summary/work/revision/epimap_overlap

module load bedtools/2.31.0

# -------------------------
# Arguments
# -------------------------
patterns="${1:?Usage: $0 <patterns> [seed]}"
SEED="${2:-42}"

# -------------------------
# Paths
# -------------------------
raw_dir="/n/scratch/users/z/ziz597/cwas/summary/data/EpiMap/epigenome_18"
SAMPLE_INFO="/n/scratch/users/z/ziz597/cwas/summary/work/revision/epimap_overlap/data/sample_info_withctDNAestimation_methods.xlsx"
GENOME_FILE="/n/scratch/users/z/ziz597/cwas/summary/data/consensus_peaks/genome.sorted"
LABEL_FILE="/n/scratch/users/z/ziz597/cwas/summary/data/EpiMap/label.csv"

H3K27AC_CF_PEAK_DIR="/n/scratch/users/z/ziz597/cwas/cfChIP_H3K27ac/data/peaks"
H3K27AC_WBC_PEAK_DIR="/n/scratch/users/z/ziz597/cwas/whole-blood/cwas/data/peaks"

out=$(echo "$patterns" | tr ',' '_')
output_dir="work/leave_one_out/${out}"
tmp_dir="${output_dir}/tmp"
mkdir -p "$output_dir" "$tmp_dir"

echo "============================================================"
echo "Leave-one-out overlap analysis"
echo "Patterns: $patterns"
echo "Output: $output_dir"
echo "============================================================"

# -------------------------
# Helpers
# -------------------------
rand_src() {
    local seed="$1"
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

# Build consensus peaks from a list of peak files
# Usage: build_consensus <file_list> <output_bed> <min_samples>
build_consensus() {
    local file_list="$1"
    local output_bed="$2"
    local min_samples="${3:-5}"

    # Cat all peaks, extend by 100bp, merge, then count overlaps
    cat $(cat "$file_list") \
        | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
        | sort -k1,1 -k2,2n \
        | bedtools slop -b 100 -g "$GENOME_FILE" -i stdin \
        | bedtools merge -i stdin \
        | sort -k1,1 -k2,2n \
        > "${output_bed}.merged.tmp"

    # Count how many samples overlap each merged region
    # For each merged region, count intersecting peak files
    local n_files=$(wc -l < "$file_list")
    bedtools intersect \
        -a "${output_bed}.merged.tmp" \
        -b $(cat "$file_list") \
        -c \
        | awk -v min="$min_samples" 'BEGIN{OFS="\t"} $4 >= min {print $1,$2,$3}' \
        > "$output_bed"

    rm -f "${output_bed}.merged.tmp"
}

# -------------------------
# Step 0: Extract sample info from xlsx using Python
# -------------------------
echo "[Step0] Extracting sample info from xlsx..."

SAMPLE_TSV="${tmp_dir}/sample_info.tsv"

python3 - "$SAMPLE_INFO" "$SAMPLE_TSV" <<'PYEOF'
import sys
import openpyxl

xlsx_file = sys.argv[1]
out_file = sys.argv[2]

wb = openpyxl.load_workbook(xlsx_file, read_only=True)
ws = wb.active

headers = [cell.value for cell in next(ws.iter_rows(min_row=1, max_row=1))]

# Find column indices
col_map = {}
for i, h in enumerate(headers):
    if h is not None:
        col_map[h.strip()] = i

needed = ['source', 'antibody', 'study_name', 'cancer_type', 'ctDNA']
for n in needed:
    if n not in col_map:
        # Try case-insensitive
        for k in col_map:
            if k.lower() == n.lower():
                col_map[n] = col_map[k]
                break

with open(out_file, 'w') as f:
    f.write("source\tantibody\tstudy_name\tcancer_type\tctDNA\n")
    for row in ws.iter_rows(min_row=2, values_only=True):
        vals = list(row)
        source = vals[col_map.get('source', 0)] or ""
        antibody = vals[col_map.get('antibody', 1)] or ""
        study_name = vals[col_map.get('study_name', 2)] or ""
        cancer_type = vals[col_map.get('cancer_type', 3)] or ""
        ctDNA = vals[col_map.get('ctDNA', 4)] or ""
        f.write(f"{source}\t{antibody}\t{study_name}\t{cancer_type}\t{ctDNA}\n")

wb.close()
PYEOF

echo "    Extracted $(tail -n +2 "$SAMPLE_TSV" | wc -l) samples"

# Filter to H3K27ac samples only
H3K27AC_SAMPLES="${tmp_dir}/h3k27ac_samples.tsv"
head -1 "$SAMPLE_TSV" > "$H3K27AC_SAMPLES"
awk -F'\t' '$2 == "H3K27ac"' "$SAMPLE_TSV" >> "$H3K27AC_SAMPLES"
echo "    H3K27ac samples: $(tail -n +2 "$H3K27AC_SAMPLES" | wc -l)"

# -------------------------
# Step 1: Build enhancer BEDs (same as original script)
# -------------------------
IFS=',' read -r -a pattern_array <<< "$patterns"
output_dir_states="${out}/1.0"
enhancer_dir="${raw_dir}/${output_dir_states}"

if ls "${enhancer_dir}"/*_18_core_K27ac_dense.bed >/dev/null 2>&1; then
    echo "[Step1] Enhancer BEDs exist, skipping build."
else
    echo "[Step1] Building enhancer BEDs..."
    mkdir -p "${enhancer_dir}"
    for file in ${raw_dir}/*.gz; do
        base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed.gz//g')
        output_file="${enhancer_dir}/H3K27ac_${base_name}_18_core_K27ac_dense.bed"
        > "$output_file"
        for pattern in "${pattern_array[@]}"; do
            zcat "$file" | grep "$pattern" | sort -k1,1 -k2,2n | bedtools merge -i - >> "$output_file"
        done
    done
fi

enh_files=($(ls "${enhancer_dir}"/*_18_core_K27ac_dense.bed))
echo "    ${#enh_files[@]} enhancer files"

# -------------------------
# Step 2: Get total counts and names for enhancer files
# -------------------------
echo "[Step2] Computing enhancer total counts and names..."

total_counts_tsv="${tmp_dir}/total_counts.tsv"
name_tsv="${tmp_dir}/name.tsv"
> "$total_counts_tsv"
> "$name_tsv"

for file in "${enh_files[@]}"; do
    base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')
    total_count=$(wc -l < "$file")
    echo -e "${base_name}\t${total_count}" >> "$total_counts_tsv"

    eid=$(echo "$base_name" | sed 's/^H3K27ac_//')
    name=$(grep -w "$eid" "$LABEL_FILE" | head -n 1 | cut -d',' -f2)
    echo -e "${base_name}\t${name}" >> "$name_tsv"
done

# -------------------------
# Step 3: Identify valid peak files for each sample
# -------------------------
echo "[Step3] Matching samples to peak files..."

CF_PEAK_LIST="${tmp_dir}/cf_peak_files.tsv"
> "$CF_PEAK_LIST"

# For each H3K27ac sample, find its peak file
while IFS=$'\t' read -r source antibody study_name cancer_type ctDNA; do
    [ "$study_name" = "study_name" ] && continue

    # Try to find peak file
    peak_file=""
    for ext in .bed .rep1_sorted_peaks.narrowPeak.bed .broadPeak; do
        candidate="${H3K27AC_CF_PEAK_DIR}/${study_name}${ext}"
        if [ -f "$candidate" ]; then
            peak_file="$candidate"
            break
        fi
    done

    if [ -n "$peak_file" ]; then
        echo -e "${study_name}\t${source}\t${cancer_type}\t${peak_file}" >> "$CF_PEAK_LIST"
    fi
done < "$H3K27AC_SAMPLES"

N_MATCHED=$(wc -l < "$CF_PEAK_LIST")
echo "    Matched $N_MATCHED cfChIP samples to peak files"

# -------------------------
# Step 4: Use pre-built WBC consensus peaks
# -------------------------
echo "[Step4] Using pre-built WBC consensus peaks..."

WBC_CONSENSUS="/n/scratch/users/z/ziz597/cwas/summary/work/revision/enrichment/data/background/wb.consensus.peaks.bed"
echo "    WBC consensus: $(wc -l < "$WBC_CONSENSUS") peaks"

# -------------------------
# Step 5: Run overlap for full dataset + leave-one-out
# -------------------------
echo "[Step5] Running overlap analyses..."

RESULTS_FILE="${output_dir}/leave_one_out_results.tsv"
echo -e "analysis\texcluded\tn_cf_samples\tn_cf_peaks\tEID\twb_count\tcf_count\ttotal_count" > "$RESULTS_FILE"

# Function to run one overlap analysis
# Usage: run_overlap <analysis_label> <excluded_label> <cf_peak_file_list>
run_overlap() {
    local label="$1"
    local excluded="$2"
    local cf_file_list="$3"
    local n_samples=$(wc -l < "$cf_file_list")

    local cf_consensus="${tmp_dir}/cf_consensus.${label}.bed"

    echo "  >>> ${label} (excluded: ${excluded}, n=${n_samples} samples)"

    # Build cfChIP consensus
    build_consensus "$cf_file_list" "$cf_consensus" 5
    local n_cf_peaks=$(wc -l < "$cf_consensus")
    echo "      cfChIP consensus: ${n_cf_peaks} peaks"

    if [ "$n_cf_peaks" -eq 0 ]; then
        echo "      [WARN] No consensus peaks, skipping."
        return
    fi

    # Remove chr prefix for consistency, filter sex/MT chroms
    local cf_nochr="${tmp_dir}/cf.nochr.${label}.bed"
    sed 's/^chr//' "$cf_consensus" | awk '$1!="X" && $1!="Y" && $1!="M"' | sort -k1,1 -k2,2n > "$cf_nochr"

    local wb_nochr="${tmp_dir}/wb.nochr.bed"
    if [ ! -f "$wb_nochr" ]; then
        sed 's/^chr//' "$WBC_CONSENSUS" | awk '$1!="X" && $1!="Y" && $1!="M"' | sort -k1,1 -k2,2n > "$wb_nochr"
    fi

    # Get per-chrom counts from cfChIP for downsampling WBC
    local target_counts="${tmp_dir}/target_counts.${label}.txt"
    cut -f1 "$cf_nochr" | sort -k1,1 | uniq -c \
        | awk 'BEGIN{OFS="\t"}{print $2,$1}' > "$target_counts"

    # Downsample WBC by chromosome
    local wb_down="${tmp_dir}/wb.down.${label}.bed"
    > "$wb_down"

    while read -r chrom cnt; do
        awk -v c="$chrom" 'BEGIN{OFS="\t"} $1==c' "$wb_nochr" \
            | shuf --random-source=<(rand_src "${SEED}_${label}_${chrom}") \
            | head -n "$cnt" \
            >> "$wb_down"
    done < "$target_counts"

    # Add chr prefix for intersect with enhancer BEDs
    local cf_chr="${tmp_dir}/cf.chr.${label}.bed"
    local wb_chr="${tmp_dir}/wb.chr.${label}.bed"
    sed 's/^/chr/' "$cf_nochr" > "$cf_chr"
    sed 's/^/chr/' "$wb_down" > "$wb_chr"

    # Overlap with each enhancer file
    for file in "${enh_files[@]}"; do
        base_name=$(basename "$file" | sed 's/_18_core_K27ac_dense.bed//g')
        total_count=$(grep -w "$base_name" "$total_counts_tsv" | cut -f2)

        wb_count=$(bedtools intersect -u -a "$file" -b "$wb_chr" | wc -l)
        cf_count=$(bedtools intersect -u -a "$file" -b "$cf_chr" | wc -l)

        echo -e "${label}\t${excluded}\t${n_samples}\t${n_cf_peaks}\t${base_name}\t${wb_count}\t${cf_count}\t${total_count}" >> "$RESULTS_FILE"
    done

    # Cleanup
    rm -f "$cf_consensus" "$cf_nochr" "$cf_chr" "$wb_down" "$wb_chr" "$target_counts"
}

# --- Full dataset ---
ALL_CF_FILES="${tmp_dir}/all_cf_files.txt"
cut -f4 "$CF_PEAK_LIST" > "$ALL_CF_FILES"
run_overlap "all" "none" "$ALL_CF_FILES"

# --- Leave-one-cancer-type-out ---
echo ""
echo ">>> Leave-one-cancer-type-out..."
cancer_types=$(cut -f3 "$CF_PEAK_LIST" | sort -u)

for ct in $cancer_types; do
    safe_ct=$(echo "$ct" | tr ' /' '_')
    ct_file_list="${tmp_dir}/cf_files.excl_ct_${safe_ct}.txt"
    awk -F'\t' -v ct="$ct" '$3 != ct {print $4}' "$CF_PEAK_LIST" > "$ct_file_list"

    n_remaining=$(wc -l < "$ct_file_list")
    if [ "$n_remaining" -lt 10 ]; then
        echo "  [SKIP] Excluding ${ct} leaves only ${n_remaining} samples"
        continue
    fi

    run_overlap "excl_cancer_${safe_ct}" "$ct" "$ct_file_list"
    rm -f "$ct_file_list"
done

# --- Leave-one-site-out ---
echo ""
echo ">>> Leave-one-site-out..."
sites=$(cut -f2 "$CF_PEAK_LIST" | sort -u)

for site in $sites; do
    safe_site=$(echo "$site" | tr ' /' '_')
    site_file_list="${tmp_dir}/cf_files.excl_site_${safe_site}.txt"
    awk -F'\t' -v s="$site" '$2 != s {print $4}' "$CF_PEAK_LIST" > "$site_file_list"

    n_remaining=$(wc -l < "$site_file_list")
    if [ "$n_remaining" -lt 10 ]; then
        echo "  [SKIP] Excluding ${site} leaves only ${n_remaining} samples"
        continue
    fi

    run_overlap "excl_site_${safe_site}" "$site" "$site_file_list"
    rm -f "$site_file_list"
done

# -------------------------
# Step 6: Summarize results in R
# -------------------------
echo ""
echo "[Step6] Summarizing results..."

Rscript - "$RESULTS_FILE" "$name_tsv" "$output_dir" <<'RSCRIPT'
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
results_file <- args[1]
name_file    <- args[2]
output_dir   <- args[3]

df <- fread(results_file)
names_df <- fread(name_file, header = FALSE, col.names = c("EID", "name"))
df <- merge(df, names_df, by = "EID", all.x = TRUE)

# Compute overlap percentages
df[, cf_pct := 100 * cf_count / total_count]
df[, wb_pct := 100 * wb_count / total_count]
df[, diff_pct := cf_pct - wb_pct]

# --- Summary table: median overlap per analysis ---
summary_dt <- df[, .(
    n_cf_samples = n_cf_samples[1],
    n_cf_peaks = n_cf_peaks[1],
    median_cf_pct = median(cf_pct),
    median_wb_pct = median(wb_pct),
    median_diff = median(diff_pct),
    n_eids = .N
), by = .(analysis, excluded)]

fwrite(summary_dt, file.path(output_dir, "leave_one_out_summary.tsv"), sep = "\t")

# --- Detailed results ---
fwrite(df, file.path(output_dir, "leave_one_out_detailed.tsv"), sep = "\t")

# --- Plot: median cfChIP-WBC difference across leave-one-out analyses ---
summary_dt$label <- ifelse(summary_dt$excluded == "none", "All samples",
                           paste0("Excl: ", summary_dt$excluded))
summary_dt$type <- ifelse(grepl("cancer", summary_dt$analysis), "Cancer type",
                   ifelse(grepl("site", summary_dt$analysis), "Site", "Full"))

p <- ggplot(summary_dt, aes(x = reorder(label, median_diff), y = median_diff, fill = type)) +
    geom_col(alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = summary_dt[analysis == "all", median_diff],
               linetype = "dashed", color = "firebrick") +
    coord_flip() +
    labs(
        title = "Median cfChIP - WBC overlap difference across EpiMap enhancers",
        subtitle = "Dashed line = full dataset",
        x = "",
        y = "Median difference in % overlap (cfChIP - WBC)",
        fill = "Exclusion type"
    ) +
    scale_fill_manual(values = c("Full" = "grey40", "Cancer type" = "#4A90D9", "Site" = "#E66100")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(output_dir, "leave_one_out_median_diff.pdf"),
       p, width = 10, height = max(6, nrow(summary_dt) * 0.4))

cat(">>> Summary saved to", file.path(output_dir, "leave_one_out_summary.tsv"), "\n")
cat(">>> Done!\n")
RSCRIPT

echo "============================================================"
echo "DONE. Results in: $output_dir"
echo "============================================================"