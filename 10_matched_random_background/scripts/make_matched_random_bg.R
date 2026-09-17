#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: make_matched_random_bg.R <fg.features.tsv> <cand.features.tsv> <out.bed>\n")
  quit(status = 1)
}

fg_path   <- args[[1]]
cand_path <- args[[2]]
out_bed   <- args[[3]]

#Test
# fg_path   <- "work/fg.features.filt.tsv"
# cand_path <- "work/cand.features.filt.tsv"
# out_bed   <- "work/out.bed"


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(nullranges)
  library(S4Vectors)
})

# ------------------------------------------------------------
# Input format (tab-separated, no header):
# chr start end name gc pctN distTSS map acc len
# ------------------------------------------------------------
read_feat <- function(path) {
  df <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr","start","end","name","gc","pctN","distTSS","map","acc","len")
  df
}

fg_df   <- read_feat(fg_path)
cand_df <- read_feat(cand_path)

# Transformations (recommended):
# - distTSS is heavy-tailed -> log1p
fg_df$logDistTSS   <- log1p(fg_df$distTSS)
cand_df$logDistTSS <- log1p(cand_df$distTSS)

# add chr as a categorical covariate (helps keep chr distribution stable)
fg_df$chrF   <- factor(fg_df$chr)
cand_df$chrF <- factor(cand_df$chr)

# Build GRanges
fg_gr <- GRanges(
  seqnames = fg_df$chr,
  ranges   = IRanges(start = fg_df$start + 1, end = fg_df$end),  # BED is 0-based start
  strand   = "*"
)
#names(fg_gr) <- fg_df$name

cand_gr <- GRanges(
  seqnames = cand_df$chr,
  ranges   = IRanges(start = cand_df$start + 1, end = cand_df$end),
  strand   = "*"
)
#names(cand_gr) <- cand_df$name

# Attach covariates
mcols(fg_gr) <- DataFrame(
  gc         = fg_df$gc,
  logDistTSS = fg_df$logDistTSS,
  map        = fg_df$map,
  acc        = fg_df$acc,
  len        = fg_df$len,
  chrF       = fg_df$chrF,
  row.names  = NULL
)

mcols(cand_gr) <- DataFrame(
  gc         = cand_df$gc,
  logDistTSS = cand_df$logDistTSS,
  map        = cand_df$map,
  acc        = cand_df$acc,
  len        = cand_df$len,
  chrF       = cand_df$chrF,
  row.names  = NULL
)
set.seed(1)

# ------------------------------------------------------------
# Matching method choice:
# - nearest neighbor without replacement is NOT supported
# - stratified sampling without replacement is the recommended alternative
# Default method is rejection sampling without replacement, but stratified is robust.
# ------------------------------------------------------------
mgr <- matchRanges(
  focal   = fg_gr,
  pool    = cand_gr,
  covar   = ~ gc + logDistTSS + map + acc + len + chrF,
  method  = "stratified",
  replace = FALSE
)

matched_gr <- matched(mgr)

# Write BED: chr start end (NO chr prefix; already nochr)
# Convert GRanges to BED 0-based half-open
options(scipen = 999)

out_df <- data.frame(
  chr   = as.character(seqnames(matched_gr)),
  start = start(matched_gr) - 1,
  end   = end(matched_gr),
  stringsAsFactors = FALSE
)

write.table(out_df, file = out_bed, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
