#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ziwei_zhang@dfci.harvard.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -t 0-3:00:00
#SBATCH -J enrichArr
#SBATCH -o logs/enrichArr_%A_%a.out
#SBATCH -e logs/enrichArr_%A_%a.err

set -euo pipefail

# -------------------------
# Inputs (passed by sbatch --export)
# -------------------------
: "${FG:?Need FG}"
: "${BG_UNIV:?Need BG_UNIV}"
: "${OUTDIR:?Need OUTDIR}"
: "${TARGET_LIST:?Need TARGET_LIST}"

# Matched-random background resources
: "${GENOME_SIZES:?Need GENOME_SIZES}"
: "${FASTA:?Need FASTA}"
: "${TSS:?Need TSS}"
: "${BLACKLIST:?Need BLACKLIST}"
: "${MAPBW:?Need MAPBW}"
: "${ACCBW:?Need ACCBW}"
: "${CAND_MULT:?Need CAND_MULT}"
: "${TMPROOT:?Need TMPROOT}"

# permutation settings
: "${NBATCHES:=50}"          # number of enrich.permute batches
: "${SEED_OFFSET:=0}"        # seed offset for reproducibility

mkdir -p logs
mkdir -p "work/output/${OUTDIR}"

# -------------------------
# Pick target for this array index
# -------------------------
TARGET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TARGET_LIST")
if [ -z "$TARGET" ]; then
  echo "[ERROR] Empty TARGET for task id ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

GROUP=$(basename "$(dirname "$TARGET")")     # eQTLs or asaQTLs
NAME=$(basename "$TARGET" .bed)             # file basename without .bed

JOBROOT="work/output/${OUTDIR}/${GROUP}/${NAME}"
mkdir -p "${JOBROOT}/fg" "${JOBROOT}/perm" "${JOBROOT}/rand"

echo "============================================================"
echo "[TASK] ${SLURM_ARRAY_TASK_ID}"
echo "[FG] $FG"
echo "[BG_UNIV] $BG_UNIV"
echo "[TARGET] $TARGET"
echo "[OUT] $JOBROOT"
echo "============================================================"

# -------------------------
# 1) Foreground enrichment
# -------------------------
sh scripts/enrich.sh "$FG" "$TARGET" "${JOBROOT}/fg/fg.enrich.txt"

# -------------------------
# 2) Background permutations (existing universe bg)
# -------------------------
for b in $(seq 1 "$NBATCHES"); do
  sh scripts/enrich.permute.sh "$b" "$BG_UNIV" "$FG" "$TARGET" "${JOBROOT}/perm"
done

# -------------------------
# 3) Matched random permutations
#    (Generate 1 matched random background per batch seed)
# -------------------------
for b in $(seq 1 "$NBATCHES"); do
  seed=$((SEED_OFFSET + b))
  bg_matched="${JOBROOT}/rand/bg.matched.seed${seed}.bed"

  scripts/make_matched_random_bg.sh \
    "$FG" \
    "$GENOME_SIZES" \
    "$FASTA" \
    "$TSS" \
    "$BLACKLIST" \
    "$MAPBW" \
    "$ACCBW" \
    "$seed" \
    "$CAND_MULT" \
    "$bg_matched" \
    "$TMPROOT"

  sh scripts/enrich.permute.sh "$b" "$bg_matched" "$FG" "$TARGET" "${JOBROOT}/rand"
done

# -------------------------
# 4) Summarize
# -------------------------
# -------------------------
# Summary + uncertainty
# -------------------------
fg_enrich=$(sed '1d' "${JOBROOT}/fg/fg.enrich.txt" | cut -f4)

# empirical p-value for BG universe null
hits=$(cat "${JOBROOT}/perm"/group.*.bg.enrich \
  | awk -v n="$fg_enrich" 'BEGIN{count=0} $1>n{count++} END{print count}')
tot=$(cat "${JOBROOT}/perm"/group.*.bg.enrich | wc -l | awk '{print $1}')

# add +1 pseudocount to avoid p=0
pval=$(echo "$hits $tot" | awk '{print ($1+1)/($2+1)}')

# mean/sd/se of null
read bg_mean bg_sd bg_se bg_n < <(cat "${JOBROOT}/perm"/group.*.bg.enrich | awk -f scripts/mean_sd_se.awk)

# empirical CI of null
bg_ci_lo=$(cat "${JOBROOT}/perm"/group.*.bg.enrich | sort -n | awk -v p=0.025 -f scripts/quantile.awk)
bg_ci_hi=$(cat "${JOBROOT}/perm"/group.*.bg.enrich | sort -n | awk -v p=0.975 -f scripts/quantile.awk)

# relative enrichment and its uncertainty
rel_enrich=$(echo "$fg_enrich $bg_mean" | awk '{print $1/$2}')

# delta-method SE for rel (fg treated fixed)
# rel = fg / bg_mean;  se(rel) ≈ rel * sd(bg)/bg_mean
rel_se=$(echo "$rel_enrich $bg_sd $bg_mean" | awk '{print $1 * ($2/$3)}')

# empirical CI for rel using rel_perm = fg / bg_perm
rel_ci_lo=$(cat "${JOBROOT}/perm"/group.*.bg.enrich \
  | awk -v fg="$fg_enrich" -v eps=1e-12 '{print fg/($1+eps)}' \
  | sort -n | awk -v p=0.025 -f scripts/quantile.awk)


rel_ci_hi=$(cat "${JOBROOT}/perm"/group.*.bg.enrich \
  | awk -v fg="$fg_enrich" -v eps=1e-12 '{print fg/($1+eps)}' \
  | sort -n | awk -v p=0.975 -f scripts/quantile.awk)

# ---- Matched-random null ----
hits_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich \
  | awk -v n="$fg_enrich" 'BEGIN{count=0} $1>n{count++} END{print count}')
tot_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich | wc -l | awk '{print $1}')
pval_rand=$(echo "$hits_rand $tot_rand" | awk '{print ($1+1)/($2+1)}')

read bg_mean_rand bg_sd_rand bg_se_rand bg_n_rand < <(cat "${JOBROOT}/rand"/group.*.bg.enrich | awk -f scripts/mean_sd_se.awk)
bg_ci_lo_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich | sort -n | awk -v p=0.025 -f scripts/quantile.awk)
bg_ci_hi_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich | sort -n | awk -v p=0.975 -f scripts/quantile.awk)

rel_enrich_rand=$(echo "$fg_enrich $bg_mean_rand" | awk '{print $1/$2}')
rel_se_rand=$(echo "$rel_enrich_rand $bg_sd_rand $bg_mean_rand" | awk '{print $1 * ($2/$3)}')

rel_ci_lo_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich \
  | awk -v fg="$fg_enrich" -v eps=1e-12 '{print fg/($1+eps)}' \
  | sort -n | awk -v p=0.025 -f scripts/quantile.awk)

rel_ci_hi_rand=$(cat "${JOBROOT}/rand"/group.*.bg.enrich \
  | awk -v fg="$fg_enrich" -v eps=1e-12 '{print fg/($1+eps)}' \
  | sort -n | awk -v p=0.975 -f scripts/quantile.awk)


# write augmented output
mkdir -p "${JOBROOT}"

printf "bg_mean\tbg_sd\tbg_se\tbg_ci_lo\tbg_ci_hi\trel_enrich\trel_se\trel_ci_lo\trel_ci_hi\tpval\tbg_mean_rand\tbg_sd_rand\tbg_se_rand\tbg_ci_lo_rand\tbg_ci_hi_rand\trel_enrich_rand\trel_se_rand\trel_ci_lo_rand\trel_ci_hi_rand\tpval_rand\n" \
  > "${JOBROOT}/enrichment.stats.txt"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$bg_mean" "$bg_sd" "$bg_se" "$bg_ci_lo" "$bg_ci_hi" \
  "$rel_enrich" "$rel_se" "$rel_ci_lo" "$rel_ci_hi" "$pval" \
  "$bg_mean_rand" "$bg_sd_rand" "$bg_se_rand" "$bg_ci_lo_rand" "$bg_ci_hi_rand" \
  "$rel_enrich_rand" "$rel_se_rand" "$rel_ci_lo_rand" "$rel_ci_hi_rand" "$pval_rand" \
  >> "${JOBROOT}/enrichment.stats.txt"

# if you still want fg columns pasted in:
paste "${JOBROOT}/fg/fg.enrich.txt" "${JOBROOT}/enrichment.stats.txt" > "${JOBROOT}/enrichment.stats.with_fg.txt"


echo "[DONE] ${GROUP}/${NAME}"
