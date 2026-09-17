#!/usr/bin/env Rscript
# Summarize Roadmap/EpiMap bootstrap overlap counts.
# Called by scripts/roadmap_overlap_bootstrapping.sh:
#   Rscript scripts/summarize_epimap_bootstrap.R \
#     <boot.long.tsv> <total_counts.tsv> <name.tsv> <combined_output.tsv>
#
# boot.long.tsv columns (no header):
#   rep  EID  WB_count  CF_count
# total_counts.tsv / name.tsv columns (no header):
#   EID  value

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: summarize_epimap_bootstrap.R <boot.long.tsv> <total_counts.tsv> <name.tsv> <combined_output.tsv>\n")
  quit(status = 1)
}

boot_path   <- args[[1]]
total_path  <- args[[2]]
name_path   <- args[[3]]
out_path    <- args[[4]]

suppressPackageStartupMessages({
  library(data.table)
})

q025 <- function(x) as.numeric(stats::quantile(x, 0.025, names = FALSE, type = 7))
q975 <- function(x) as.numeric(stats::quantile(x, 0.975, names = FALSE, type = 7))

boot <- fread(
  boot_path,
  header = FALSE,
  col.names = c("rep", "EID", "WB_count", "CF_count")
)
totals <- fread(total_path, header = FALSE, col.names = c("EID", "Total_count"))
names_df <- fread(name_path, header = FALSE, col.names = c("EID", "NAME"))

summary_dt <- boot[, .(
  n_boot     = .N,
  WB_count   = mean(WB_count),
  WB_sd      = sd(WB_count),
  WB_ci_lo   = q025(WB_count),
  WB_ci_hi   = q975(WB_count),
  CF_count   = mean(CF_count),
  CF_sd      = sd(CF_count),
  CF_ci_lo   = q025(CF_count),
  CF_ci_hi   = q975(CF_count)
), by = EID]

summary_dt <- merge(summary_dt, totals, by = "EID", all.x = TRUE)
summary_dt <- merge(summary_dt, names_df, by = "EID", all.x = TRUE)

# Prefix matches 01_data_preprocessing combined_results.tsv
# (EID, WB_count, CF_count, Total_count, NAME) so existing plots can reuse means.
out <- summary_dt[, .(
  EID, WB_count, CF_count, Total_count, NAME,
  n_boot, WB_sd, WB_ci_lo, WB_ci_hi, CF_sd, CF_ci_lo, CF_ci_hi
)]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_path, sep = "\t")

cat("Wrote bootstrap summary:", out_path, "\n")
cat("  n epigenomes:", nrow(out), "\n")
