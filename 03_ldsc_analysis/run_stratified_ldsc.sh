#!/bin/bash
#SBATCH --job-name=s-ldsc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ziwei_zhang@dfci.harvard.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -t 0-3



trait=$1
models=$2
baseline_model=$3

cd /n/scratch/users/z/ziz597/cwas/summary/work/S-ldsc
mkdir -p work/heritability_enrichment/${trait}
python ldsc/ldsc.py --h2 work/sumstats/${trait}.sumstats.gz --w-ld-chr data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --ref-ld-chr $models --overlap-annot --frqfile-chr data/1000G_Phase3_frq/1000G.EUR.QC. --out work/heritability_enrichment/${trait}/${trait}_${baseline_model} --print-coefficients



