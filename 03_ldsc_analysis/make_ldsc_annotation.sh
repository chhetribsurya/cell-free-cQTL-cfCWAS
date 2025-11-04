#!/bin/bash
#SBATCH --job-name=make_annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ziwei_zhang@dfci.harvard.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -t 0-3



bed=$1
name=$2
out=$3

cd /n/scratch/users/z/ziz597/cwas/summary/work/S-ldsc

mkdir -p work/models/${out}
sort -k1,1n -k2,2n $bed | sed 's/^/chr/g' | cut -f 1-3 | sed 's/ /\t/g' > work/models/${out}/out.bed
for chr in {1..22}
do
	python ldsc/make_annot.py --bed-file work/models/${out}/out.bed --bimfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim --annot-file work/models/${out}/${name}.${chr}.annot.gz
	python ldsc/ldsc.py --l2 --bfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} --ld-wind-cm 1 --annot work/models/${out}/${name}.${chr}.annot.gz --thin-annot --out work/models/${out}/${name}.${chr} --print-snps data/hapmap3_snps/list.txt 
done

