# ---- Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))
setwd("~/Projects/cfChIP")

# ---- Inputs
peaks <- c("k27_combined_cQTLs_Gusev_asqtls",
           "k4_combined_cQTLs_Gusev_asqtls",
           "wb_combined_cQTLs_Gusev_asqtls")

# ---- Load & tidy
dat <- lapply(peaks, function(peak) {
  fn <- file.path("Results/Gusev_asQTLs_enrichment", paste0(peak, ".enrichment.summary.txt"))
  d  <- fread(fn, sep = "\t", header = FALSE, colClasses = c("character", rep("numeric", 10)))
  setnames(d, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11"),
           c("Tissues","fg_tot","fg_hit","fg_num","fg_enrich","bg_enrich",
             "rel_enrich","pval","bg_enrich_rand","rel_enrich_rand","pval_rand"))
  d$annot <- peak
  d
}) |> rbindlist()

# Friendly labels & fixed order
dat <- dat |>
  mutate(
    annot = factor(annot, levels = peaks,
                   labels = c("H3K27ac \ncfcQTLs", "H3K4me3 \ncfcQTLs", "H3K27ac \nWBC cQTLs"))
  )

# ---- Significance filter (Bonferroni for 0.1 across 72 tests)
ALPHA <- 0.05/72
boxplot.df <- dat |> filter(pval < ALPHA)

# Optional quick counts
# table(boxplot.df$annot)

# ---- Plot
p <- ggplot(boxplot.df, aes(x = annot, y = rel_enrich)) +
  geom_boxplot(fill = "#1f78b4", color = "black", outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.25, alpha = 0.7, size = 1.5) +
  labs(x = NULL,
       y = "Fold enrichment at cancer as-aQTLs") +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    plot.margin  = margin(5, 10, 5, 5)
  )

print(p)

# ---- Save vector PDF (no dpi needed)
ggsave(
  filename = "Manuscripts/Figures/3_f.pdf",
  plot = p,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  width = 100, height = 140, units = "mm"
)

# ---- Supplement table: wide format
dat_wide <- dat |>
  select(Tissues, rel_enrich, annot) |>
  pivot_wider(names_from = annot, values_from = rel_enrich)

write_excel_csv(dat_wide, "Tables/sup8.csv")
