# ---- Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))
setwd("~/Projects/cfChIP")

# ---- Load
traits <- c(
  "Mean_platelet_volume",
  "Plateletcrit",
  "UKB_460K.blood_PLATELET_COUNT",
  "UKB_460K.blood_PLATELET_DISTRIB_WIDTH"
)

dat <- lapply(traits, function(trait) {
  fn <- file.path("Results/S_LDSC", paste0(trait, ".ldsc.enrichment.txt"))
  # Expect columns: Peaks, (skip 3), enrichment, enrichment.se, enrichment.p-value, (skip 3)
  d  <- fread(fn, sep = "\t", header = FALSE,
              colClasses = c("character", rep("NULL", 3), rep("numeric", 3), rep("NULL", 3)))
  setnames(d, c("V1","V5","V6","V7"), c("Peaks","enrichment","enrichment.se","enrichment.p"))
  d$trait <- trait
  d
}) |>
  rbindlist()

# ---- Keep the three peak sets and relabel
dat <- dat |>
  filter(Peaks %in% c("cfChIP_H3K27ac.combined.sig",
                      "cfChIP_H3K4me3.combined.sig.peak",
                      "whole-blood.combined.sig.peak")) |>
  mutate(
    Peaks = case_when(
      Peaks == "cfChIP_H3K27ac.combined.sig"        ~ "H3K27ac \ncfcQTLs",
      Peaks == "cfChIP_H3K4me3.combined.sig.peak"   ~ "H3K4me3 \ncfcQTLs",
      Peaks == "whole-blood.combined.sig.peak"      ~ "H3K27ac \nWBC cQTLs",
      TRUE ~ Peaks
    ),
    trait = case_when(
      trait == "Mean_platelet_volume"                   ~ "MPV",
      trait == "Plateletcrit"                           ~ "PCT",
      trait == "UKB_460K.blood_PLATELET_COUNT"          ~ "PLT",
      trait == "UKB_460K.blood_PLATELET_DISTRIB_WIDTH"  ~ "PDW",
      TRUE ~ trait
    )
  )

# ---- Significance stars
dat <- dat |>
  mutate(sig = case_when(
    enrichment.p <= 0.001 ~ "***",
    enrichment.p <= 0.01  ~ "**",
    enrichment.p <= 0.05  ~ "*",
    TRUE                  ~ "ns"
  ))

# ---- Ordering (same spirit as your script)
# Order peaks by enrichment within MPV (desc)
peak_levels <- dat |> filter(trait == "MPV") |>
  arrange(desc(enrichment)) |> pull(Peaks) |> unique()
dat$Peaks <- factor(dat$Peaks, levels = peak_levels)

# Order traits by enrichment within H3K27ac cf-cQTLs (desc)
trait_levels <- dat |> filter(Peaks == "H3K27ac \ncfcQTLs") |>
  arrange(desc(enrichment)) |> pull(trait) |> unique()
dat$trait <- factor(dat$trait, levels = trait_levels)

# ---- Plot
dodge_w <- 0.6
p <- ggplot(dat, aes(x = trait, y = enrichment, fill = Peaks, group = Peaks)) +
  geom_col(position = position_dodge(width = dodge_w), width = dodge_w, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = enrichment - enrichment.se, ymax = enrichment + enrichment.se),
                width = 0.12, position = position_dodge(width = dodge_w)) +
  geom_text(aes(label = sig),
            position = position_dodge(width = dodge_w),
            vjust = -5, hjust = -0.1, size = 10/.pt, colour = "black") +
  scale_fill_manual(
    values = c("H3K4me3 \ncfcQTLs" = "#50AD9F",
               "H3K27ac \ncfcQTLs" = "#E9C716",
               "H3K27ac \nWBC cQTLs"        = "#B22222"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Platelet traits",
    y = "Heritability enrichment"
  ) +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 10, colour = "black"),
    axis.text.x  = element_text(size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    plot.margin  = margin(5, 15, 5, 5)
  )

print(p)

# ---- Save vector PDF (no DPI needed)
ggsave(
  filename = "Manuscripts/Figures/3_h.pdf",
  plot = p,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  width = 100, height = 140, units = "mm"
)
