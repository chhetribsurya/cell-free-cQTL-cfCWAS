#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your.email@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -t 0-1

fg=$1
bg=$2
line=$3
target=$4
out_dir=$5

#fg
mkdir -p work/enrichment/${out_dir}/$line
sh scripts/enrich.sh $fg $target work/enrichment/${out_dir}/${line}/fg.enrich.txt

#permute
mkdir -p work/enrichment/${out_dir}/perm/$line
for i in $(seq 1 50)
do 
	sh scripts/enrich.permute.sh $i $bg $fg $target work/enrichment/${out_dir}/perm/$line
done
#
##random
mkdir -p work/enrichment/${out_dir}/perm/${line}/rand
for i in $(seq 1 50)
do
	printf "" > work/enrichment/${out_dir}/perm/${line}/rand/tmp.${i}
	bedtools shuffle -seed $i -chrom -noOverlapping -i $fg -g /path/to/reference/hg19/hg19.chrom.sizes.nochr  >> work/enrichment/${out_dir}/perm/${line}/rand/tmp.${i}
	sh scripts/enrich.permute.sh ${i} work/enrichment/${out_dir}/perm/${line}/rand/tmp.${i} $fg $target work/enrichment/${out_dir}/perm/${line}/rand
	rm work/enrichment/${out_dir}/perm/${line}/rand/tmp.${i}
done

#Summary
#mkdir work/enrichment/results
fg_enrich=`sed '1d' work/enrichment/${out_dir}/${line}/fg.enrich.txt | cut -f4`
hits=`cat work/enrichment/${out_dir}/perm/${line}/group.*.bg.enrich | awk -v n=$fg_enrich 'BEGIN{{count=0}} $1>n{{count++}}END{{print count}}'`
tot=`cat work/enrichment/${out_dir}/perm/${line}/group.*.bg.enrich | wc -l`
bg_enrich=`cat work/enrichment/${out_dir}/perm/${line}/group.*.bg.enrich | awk -f scripts/avg.awk`
rel_enrich=`echo $fg_enrich $bg_enrich | awk '{{print $1/$2}}'`
pval=`echo $hits $tot | awk '{{print $1/$2}}'`
hits_rand=`cat work/enrichment/${out_dir}/perm/${line}/rand/group.*.bg.enrich | awk -v n=$fg_enrich 'BEGIN{{count=0}} $1>n{{count++}}END{{print count}}'`
tot_rand=`cat work/enrichment/${out_dir}/perm/${line}/rand/group.*.bg.enrich | wc -l`
bg_enrich_rand=`cat work/enrichment/${out_dir}/perm/${line}/rand/group.*.bg.enrich | awk -f scripts/avg.awk`
rel_enrich_rand=`echo $fg_enrich $bg_enrich_rand | awk '{{print $1/$2}}'`
pval_rand=`echo $hits_rand $tot_rand | awk '{{print $1/$2}}'`
mkdir -p work/enrichment/${out_dir}/results/${line}
printf "bg_enrich\trel_enrich\tpval\tbg_enrich_rand\trel_enrich_rand\tpval_rand\n%s\t%s\t%s\t%s\t%s\t%s\n" $bg_enrich $rel_enrich $pval $bg_enrich_rand $rel_enrich_rand $pval_rand| paste work/enrichment/${out_dir}/${line}/fg.enrich.txt - > work/enrichment/${out_dir}/results/${line}/enrichment.stats.txt

#plot
#grep -v fg work/enrichment/results/*/enrichment.stats.txt | sed 's/.*results\///' | sed 's/\/enrichment.stats.txt:/\t/' > work/enrichment/results/enrichment.summary.txt
#Rscript scripts/plot_enrichment.R work/enrichment/results/enrichment.summary.txt work/enrichment/results/enrichment.pdf
