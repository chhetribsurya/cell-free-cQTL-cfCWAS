#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# make_matched_random_bg.sh
#
# Generate genome-wide random background matched to FG on:
#   - chr, length (by construction)
#   - GC, distance-to-TSS, mappability, accessibility (by matchRanges)
#
# Output: bg.matched.bed  (no header, 3 columns, NO 'chr' prefix)
# ============================================================

if [ $# -lt 10 ]; then
  echo "Usage:"
  echo "  $0 <fg.bed> <genome.sizes.nochr> <fasta.fa> <tss.bed.nochr> <blacklist.bed.nochr> <mappability.bw> <accessibility.bw> <seed> <multiplier> <out.bed>"
  echo ""
  echo "Example:"
  echo "  $0 fg.bed /home/.../hg19.chrom.sizes.nochr /home/.../hg19.fa tss.nochr.bed blacklist.nochr.bed map.bw acc.bw 1 200 bg.matched.bed"
  exit 1
fi

FG=$1
GENOME_SIZES=$2
FASTA=$3
TSS_BED=$4
BLACKLIST=$5
MAP_BW=$6
ACC_BW=$7
SEED=$8
MULT=$9
OUT_BED=${10}

# temp workspace
WORKDIR=$(mktemp -d -p "${TMPDIR:-/tmp}" matchedbg.XXXXXX)
trap 'rm -rf "$WORKDIR"' EXIT

# ---------------------------
# 1) Standardize FG (no chr prefix, add unique ID)
# ---------------------------
FG_ID=${WORKDIR}/fg.id.bed
cat "$FG" \
  | sed 's/^chr//' \
  | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"fg_"NR}' \
  > "$FG_ID"

FG_N=$(wc -l < "$FG_ID")
CAND_N=$((FG_N * MULT))

# ---------------------------
# 2) Make candidate pool (chr + length matched)
#    by replicating FG then shuffling genome-wide
# ---------------------------
FG_REP=${WORKDIR}/fg.rep.bed
CAND_ID=${WORKDIR}/cand.id.bed

awk -v m="$MULT" 'BEGIN{OFS="\t"}{for(i=1;i<=m;i++) print $1,$2,$3,$4"_"i}' "$FG_ID" > "$FG_REP"

bedtools shuffle \
  -i "$FG_REP" \
  -g "$GENOME_SIZES" \
  -chrom \
  -excl "$BLACKLIST" \
  -seed "$SEED" \
  | sort -k1,1 -k2,2n > "$CAND_ID"

# ---------------------------
# 3) Prepare sorted TSS
# ---------------------------
TSS_SORT=${WORKDIR}/tss.sorted.bed
cat "$TSS_BED" | sed 's/^chr//' | sort -k1,1 -k2,2n > "$TSS_SORT"

# ============================================================
# Helper: extract named columns from bedtools nuc
# ============================================================
extract_nuc_cols () {
  local nuc_in=$1
  local out_tsv=$2
  # output: name \t pct_gc \t pct_n
  awk 'BEGIN{OFS="\t"} NR>1{
  id=$4
  pct_gc=$6
  pct_n=($13>0 ? $11/$13 : 0)
  print id, pct_gc, pct_n
  }' "$nuc_in" > "$out_tsv"
}

# ---------------------------
# 4) Compute covariates: GC + pctN
# ---------------------------
FG_NUC=${WORKDIR}/fg.nuc.txt
CAND_NUC=${WORKDIR}/cand.nuc.txt
FG_GC=${WORKDIR}/fg.gc.tsv
CAND_GC=${WORKDIR}/cand.gc.tsv

FG_ID_CHR=${WORKDIR}/fg.id.chr.bed
CAND_ID_CHR=${WORKDIR}/cand.id.chr.bed

awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4}' "$FG_ID"   > "$FG_ID_CHR"
awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4}' "$CAND_ID" > "$CAND_ID_CHR"

bedtools nuc -fi "$FASTA" -bed "$FG_ID_CHR"   > "$FG_NUC"
bedtools nuc -fi "$FASTA" -bed "$CAND_ID_CHR" > "$CAND_NUC"

extract_nuc_cols "$FG_NUC"   "$FG_GC"
extract_nuc_cols "$CAND_NUC" "$CAND_GC"

# ---------------------------
# 5) Distance to nearest TSS
# ---------------------------
FG_DIST=${WORKDIR}/fg.dist.tsv
CAND_DIST=${WORKDIR}/cand.dist.tsv

bedtools closest -a "$FG_ID"   -b "$TSS_SORT" -d \
  | awk 'BEGIN{OFS="\t"}{print $4,$NF}' > "$FG_DIST"

bedtools closest -a "$CAND_ID" -b "$TSS_SORT" -d \
  | awk 'BEGIN{OFS="\t"}{print $4,$NF}' > "$CAND_DIST"

# ---------------------------
# 6) Mappability + Accessibility via bigWigAverageOverBed
#    Use mean0 (column 5): includes uncovered bases as zero.
# ---------------------------
FG_MAPTAB=${WORKDIR}/fg.map.tab
CAND_MAPTAB=${WORKDIR}/cand.map.tab
FG_ACCTAB=${WORKDIR}/fg.acc.tab
CAND_ACCTAB=${WORKDIR}/cand.acc.tab

FG_MAP=${WORKDIR}/fg.map.tsv
CAND_MAP=${WORKDIR}/cand.map.tsv
FG_ACC=${WORKDIR}/fg.acc.tsv
CAND_ACC=${WORKDIR}/cand.acc.tsv

bigWigAverageOverBed "$MAP_BW" "$FG_ID_CHR"   "$FG_MAPTAB"
bigWigAverageOverBed "$MAP_BW" "$CAND_ID_CHR" "$CAND_MAPTAB"
bigWigAverageOverBed "$ACC_BW" "$FG_ID_CHR"   "$FG_ACCTAB"
bigWigAverageOverBed "$ACC_BW" "$CAND_ID_CHR" "$CAND_ACCTAB"

# name \t mean0
awk 'BEGIN{OFS="\t"}{print $1,$5}' "$FG_MAPTAB"   > "$FG_MAP"
awk 'BEGIN{OFS="\t"}{print $1,$5}' "$CAND_MAPTAB" > "$CAND_MAP"
awk 'BEGIN{OFS="\t"}{print $1,$5}' "$FG_ACCTAB"   > "$FG_ACC"
awk 'BEGIN{OFS="\t"}{print $1,$5}' "$CAND_ACCTAB" > "$CAND_ACC"

# ---------------------------
# 7) Length
# ---------------------------
FG_LEN=${WORKDIR}/fg.len.tsv
CAND_LEN=${WORKDIR}/cand.len.tsv

awk 'BEGIN{OFS="\t"}{print $4,($3-$2)}' "$FG_ID"   > "$FG_LEN"
awk 'BEGIN{OFS="\t"}{print $4,($3-$2)}' "$CAND_ID" > "$CAND_LEN"

# ---------------------------
# 8) Assemble feature tables (by ID join)
#    Final format:
#      chr start end name gc pctN distTSS map acc len
# ---------------------------
FG_FEAT=${WORKDIR}/fg.features.tsv
CAND_FEAT=${WORKDIR}/cand.features.tsv

# Start with BED coords
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' "$FG_ID"   > ${WORKDIR}/fg.coord.tsv
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' "$CAND_ID" > ${WORKDIR}/cand.coord.tsv

# join helper: join on "name" (4th col in coord file vs 1st col in covariate files)
# We'll reorder to: name first for joining, then restore.
reorder_name_first () {
  local in_tsv=$1
  local out_tsv=$2
  # chr start end name  -> name chr start end
  awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3}' "$in_tsv" > "$out_tsv"
}

reorder_name_first ${WORKDIR}/fg.coord.tsv   ${WORKDIR}/fg.coord.namefirst.tsv
reorder_name_first ${WORKDIR}/cand.coord.tsv ${WORKDIR}/cand.coord.namefirst.tsv

# sort all join inputs
sort -k1,1 ${WORKDIR}/fg.coord.namefirst.tsv   > ${WORKDIR}/fg.coord.sorted.tsv
sort -k1,1 ${WORKDIR}/cand.coord.namefirst.tsv > ${WORKDIR}/cand.coord.sorted.tsv
sort -k1,1 "$FG_GC"   > ${WORKDIR}/fg.gc.sorted.tsv
sort -k1,1 "$CAND_GC" > ${WORKDIR}/cand.gc.sorted.tsv
sort -k1,1 "$FG_DIST"   > ${WORKDIR}/fg.dist.sorted.tsv
sort -k1,1 "$CAND_DIST" > ${WORKDIR}/cand.dist.sorted.tsv
sort -k1,1 "$FG_MAP"   > ${WORKDIR}/fg.map.sorted.tsv
sort -k1,1 "$CAND_MAP" > ${WORKDIR}/cand.map.sorted.tsv
sort -k1,1 "$FG_ACC"   > ${WORKDIR}/fg.acc.sorted.tsv
sort -k1,1 "$CAND_ACC" > ${WORKDIR}/cand.acc.sorted.tsv
sort -k1,1 "$FG_LEN"   > ${WORKDIR}/fg.len.sorted.tsv
sort -k1,1 "$CAND_LEN" > ${WORKDIR}/cand.len.sorted.tsv

# FG join chain
join -t $'\t' -1 1 -2 1 ${WORKDIR}/fg.coord.sorted.tsv ${WORKDIR}/fg.gc.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/fg.dist.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/fg.map.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/fg.acc.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/fg.len.sorted.tsv \
  | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10}' \
  > "$FG_FEAT"

# Candidate join chain
join -t $'\t' -1 1 -2 1 ${WORKDIR}/cand.coord.sorted.tsv ${WORKDIR}/cand.gc.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/cand.dist.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/cand.map.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/cand.acc.sorted.tsv \
  | join -t $'\t' -1 1 -2 1 - ${WORKDIR}/cand.len.sorted.tsv \
  | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10}' \
  > "$CAND_FEAT"

# ---------------------------
# 9) Filter out poor candidates (high N fraction)
# ---------------------------
FG_FEAT_FILT=${WORKDIR}/fg.features.filt.tsv
CAND_FEAT_FILT=${WORKDIR}/cand.features.filt.tsv

# cols: chr start end name gc pctN distTSS map acc len
awk 'BEGIN{OFS="\t"} $6 < 0.05 {print}' "$FG_FEAT" > "$FG_FEAT_FILT"
awk 'BEGIN{OFS="\t"} $6 < 0.05 {print}' "$CAND_FEAT" > "$CAND_FEAT_FILT"

# ---------------------------
# 10) Run matchRanges in R and write OUT_BED
# ---------------------------
Rscript scripts/make_matched_random_bg.R \
  "$FG_FEAT_FILT" \
  "$CAND_FEAT_FILT" \
  "$OUT_BED"

# final sort/uniq for safety
sort -k1,1 -k2,2n "$OUT_BED" | uniq > "${OUT_BED}.tmp"
mv "${OUT_BED}.tmp" "$OUT_BED"

echo "[OK] Wrote matched background: $OUT_BED"
echo "     (seed=$SEED, MULT=$MULT, fg_n=$FG_N, cand_n=$CAND_N)"
