# ---- Libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/Projects/CWAS/cfChIP")

# ---- Input
summ_file <- "Results/MK_cells_enrichment/MK_eQTLs_enrichment.csv"

# Expect columns including: Assays, Cell_type, rel_enrich, pval, rel_enrich_rand, pval_rand
d <- read_csv(summ_file, show_col_types = FALSE) %>%
  select(Assays, Cell_type, rel_enrich, pval, rel_enrich_rand, pval_rand)

# Friendly labels for assays (adjust if your raw names differ)
assay_map <- c(
  "cfChIP_K4"  = "H3K4me3 \ncf-cQTLs",
  "cfChIP_K27" = "H3K27ac \ncf-cQTLs",
  "wb_ChIP_K27"  = "H3K27ac \nWBC cQTLs"
)
d <- d %>%
  mutate(Assays = recode(Assays, !!!assay_map, .default = Assays),
         Cell_type = factor(Cell_type, levels = c("Platelets","Megakaryocytes")))

# We'll plot the *randomized background* enrichment with its p-value
rand_df <- d %>%
  transmute(Assays, Cell_type, value = rel_enrich_rand, p = pval_rand)

# Order assays by Platelets enrichment (desc); fallback to overall mean if Platelets missing
assay_levels <- rand_df %>%
  filter(Cell_type == "Platelets") %>%
  arrange(desc(value)) %>%
  pull(Assays) %>%
  unique()

if (length(assay_levels) == 0) {
  assay_levels <- rand_df %>%
    group_by(Assays) %>%
    summarise(m = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Assays)
}
rand_df <- rand_df %>% mutate(Assays = factor(Assays, levels = assay_levels))

# ---- Plotmath label helper for p-values (vectorized)
# Returns strings like "p == 1.23 %*% 10^-4" or "p < 2.2 %*% 10^-16" or "p == 0.07"
# ---- Modified p-value label function for "***"
p_to_plotmath_vec <- function(p, thresh = .Machine$double.eps) {
  out <- rep(NA_character_, length(p))
  na_mask <- is.na(p)
  
  # "***" for p < 0.001
  stars <- !na_mask & p < 0.001
  out[stars] <- "***"
  
  # Plain decimals for p > 0.05
  plain <- !na_mask & p > 0.05
  out[plain] <- paste0("p == ", sprintf("%.2f", p[plain]))
  
  # Very small p-values
  tiny <- !na_mask & !plain & !stars & (p < thresh | p == 0)
  out[tiny] <- "p < 2.2 %*% 10^-16"
  
  # Scientific notation for the rest
  sci <- !na_mask & !plain & !stars & !tiny
  if (any(sci)) {
    e <- floor(log10(p[sci]))
    m <- p[sci] / (10^e)
    out[sci] <- paste0("p == ", sprintf("%.2f", m), " %*% 10^", e)
  }
  out
}

rand_df <- rand_df %>%
  mutate(p_label = p_to_plotmath_vec(p))

# Position labels above bars
label_pos <- rand_df %>%
  group_by(Assays, Cell_type) %>%
  summarise(y = max(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y + 0.03 * max(y, na.rm = TRUE))

ann <- rand_df %>%
  left_join(label_pos, by = c("Assays","Cell_type")) %>%
  select(Assays, Cell_type, p_label, y)

# ---- Plot
p <- ggplot(rand_df, aes(x = Assays, y = value, fill = Cell_type)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6, colour = "black", linewidth = 0.2) +
  geom_text(
    data = ann,
    aes(x = Assays, y = y, label = p_label, group = Cell_type),
    position = position_dodge(width = 0.6),
    vjust = 0, size = 10/.pt, colour = "black", parse = FALSE
  ) +
  scale_fill_manual(values = c("Megakaryocytes" = "#1f78b4", "Platelets" = "#e31a1c"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Fold enrichment at eQTLs"
  ) +
  theme(
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    legend.text  = element_text(size = 10, colour = "black"),
    legend.position = "bottom"
  )

print(p)

# ---- Save vector PDF (journal-friendly)
ggsave(
  filename = "Manuscripts/Figures/3_g.pdf",
  plot = p,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  width = 100, height = 140, units = "mm"
)
