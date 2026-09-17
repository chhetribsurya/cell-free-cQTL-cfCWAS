conda activate enrich-matched-bg

FG=data/foreground/cfChIP_H3K4me3.combined.sig.bed
BG_UNIV=data/background/cfChIP_H3K4me3.consensus.peaks.bed
TARGET_DIR=data/target

OUTDIR=run_cfChIP_H3K4me3
TARGET_LIST=work/target_list.txt

GENOME_SIZES=/home/ziz597/ref/hg19/hg19.chrom.sizes.nochr
FASTA=/home/ziz597/ref/hg19/hg19.fa

TSS=data/resources/refGene_hg19_TSS.bed
BLACKLIST=data/resources/hg19-blacklist.v2.bed

MAPBW=data/resources/wgEncodeDukeMapabilityUniqueness35bp.bigWig
ACCBW=data/resources/wgEncodeCrgMapabilityAlign100mer.bigWig


NBATCHES=50          # number of group.* batches
CAND_MULT=50         # candidate multiplier used inside matched random bg
TMPROOT=/tmp         # fast local scratch
SEED_OFFSET=0

find $TARGET_DIR -type f -name "*.bed" | sort > $TARGET_LIST
wc -l $TARGET_LIST
head -n 5 $TARGET_LIST

NTARGETS=$(wc -l < "$TARGET_LIST")
echo "$NTARGETS targets"

sbatch --array=1 \
  --export=ALL,FG="$FG",BG_UNIV="$BG_UNIV",OUTDIR="$OUTDIR",TARGET_LIST="$TARGET_LIST",NBATCHES="$NBATCHES",SEED_OFFSET="$SEED_OFFSET",GENOME_SIZES="$GENOME_SIZES",FASTA="$FASTA",TSS="$TSS",BLACKLIST="$BLACKLIST",MAPBW="$MAPBW",ACCBW="$ACCBW",CAND_MULT="$CAND_MULT",TMPROOT="$TMPROOT" \
  scripts/enrichment_array.sh