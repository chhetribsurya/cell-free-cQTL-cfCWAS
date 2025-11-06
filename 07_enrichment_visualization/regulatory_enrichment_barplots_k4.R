# --- Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))
setwd("~/Projects/cfChIP/Manuscripts/Figures/")

# --- Inputs
fusion.path <- '../../Results/CWAS/Results/cfChIP_H3K4me3/Non-Cancer/'
trait.list  <- list.files(fusion.path)

# Dummy SNPs to preserve full chromosome spans (optional but helpful)
dummies <- readRDS("../../Results/dummysnps.rds")  # must have SNP, CHR, BP, P columns

# --- Collect data
plot.dat <- NULL
sig_tbl  <- list()

for (trait in trait.list) {
  cv_file  <- file.path(fusion.path, trait, "cv",  "cwas.cv.best.txt")
  sig_file <- file.path(fusion.path, trait, "sig", "cwas.sig.best.txt")
  if (!file.exists(cv_file)) next
  
  model <- fread(cv_file, sep = "\t", header = TRUE)
  if (!all(c("BEST.GWAS.ID","CHR","P0","P1","TWAS.P","TWAS.Z","BEST.GWAS.Z") %in% names(model))) next
  
  sigline <- 0
  if (file.exists(sig_file)) {
    model.sig <- tryCatch(fread(sig_file, sep = "\t", header = TRUE), error = function(e) data.table())
    if (nrow(model.sig) > 0 && "TWAS.P" %in% names(model.sig)) {
      sigline <- max(model.sig$TWAS.P, na.rm = TRUE)
      if (!is.finite(sigline)) sigline <- 0
    }
  }
  sig_tbl[[trait]] <- sigline
  
  df <- data.frame(
    SNP    = model$BEST.GWAS.ID,
    CHR    = as.integer(model$CHR),
    BP     = model$P0 + (model$P1 - model$P0)/2,
    P      = as.numeric(model$TWAS.P),
    zscore = as.numeric(model$TWAS.Z),
    gwas.z = as.numeric(model$BEST.GWAS.Z),
    GWAS   = trait,
    stringsAsFactors = FALSE
  )
  
  # pad with dummies for this trait
  d.tmp <- dummies |>
    transmute(
      SNP, CHR = as.integer(CHR), BP, P,
      zscore = NA_real_, gwas.z = NA_real_, GWAS = trait
    )
  plot.dat <- bind_rows(plot.dat, df, d.tmp)
}

stopifnot(nrow(plot.dat) > 0)

# --- Order traits by #signals (ascending)
plot.dat <- plot.dat |>
  group_by(GWAS) |>
  mutate(
    sigline = sig_tbl[[unique(GWAS)]],
    sig     = ifelse(is.finite(P), P <= sigline, NA)
  ) |>
  ungroup()%>%filter(GWAS!="Psoriasis")

trait_order <- plot.dat |>
  group_by(GWAS) |>
  summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") |>
  arrange(n_sig) |>
  pull(GWAS)

plot.dat$GWAS <- factor(plot.dat$GWAS, levels = trait_order)

# --- Compute cumulative genomic position for Manhattan
# Get per-chromosome max BP to build cumulative offsets
chr_info <- plot.dat |>
  filter(!is.na(CHR), !is.na(BP)) |>
  group_by(CHR) |>
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") |>
  arrange(CHR) |>
  mutate(add = lag(cumsum(chr_len), default = 0))

plot.dat <- plot.dat |>
  left_join(chr_info, by = "CHR") |>
  mutate(pos_cum = BP + add)

# X-axis labels at chromosome centers
axis_df <- chr_info |>
  mutate(center = add + chr_len / 2)
axis_df$label <- ifelse(axis_df$CHR %% 2 == 0, axis_df$CHR, "")

# Manhattan colors
base_col   <- "grey70"
sig_cols   <- brewer.pal(max(3, min(8, nlevels(plot.dat$GWAS))), "Dark2")
names(sig_cols) <- levels(plot.dat$GWAS)

# Per-trait horizontal threshold dataframe
thr_df <- plot.dat |>
  distinct(GWAS, sigline) |>
  filter(is.finite(sigline)) |>
  mutate(y = -log10(sigline))
thr_df = thr_df[order(thr_df$y), ][1,]

# --- Manhattan (ggplot)
# Colors per trait
trait_cols <- RColorBrewer::brewer.pal(max(3, min(8, nlevels(plot.dat$GWAS))), "Dark2")
names(trait_cols) <- levels(plot.dat$GWAS)

# One-panel Manhattan: all traits overlap
# helpers
pts_ns  <- dplyr::filter(plot.dat, is.finite(P), !sig)   # non-sig
pts_sig <- dplyr::filter(plot.dat, is.finite(P),  sig)   # sig

manh_overlap <- ggplot(plot.dat, aes(x = pos_cum)) +
  # non-significant points (grey, no legend)
  geom_point(data = pts_ns,
             aes(y = -log10(P)),
             colour = "grey75", size = 1, alpha = 0.5, show.legend = FALSE) +
  # significant points (colored by trait/panel, keep legend)
  geom_point(data = pts_sig,
             aes(y = -log10(P), colour = GWAS),
             size = 1, alpha = 0.8) +
  # per-trait threshold lines
  geom_hline(data = thr_df,
             aes(yintercept = -log10(sigline)), colour = 'red',
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_color_manual(values = trait_cols, name = NULL) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$label,
                     expand = expansion(mult = c(0.002, 0.01))) +
  scale_y_continuous(name = expression(-log[10](P)),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", title = "H3K4me3 cfCWAS") +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 10, colour = "black"),
    axis.text    = element_text(size = 8,  colour = "black"),
    axis.title   = element_text(size = 10, colour = "black"),
    plot.margin  = margin(5, 20, 5, 5)
  )

manh_overlap

ggsave("4_b.pdf", manh_overlap, device = "pdf",
       width = 140/25.4, height = 100/25.4, units = "in")
