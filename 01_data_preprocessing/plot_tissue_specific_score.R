setwd("~/Projects/cfChIP")
# visualize_detection.R
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)   # for percent()
library(RColorBrewer)
library(extrafont)  # For managing fonts
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 7)
)

# 1) Read in the counts
bg.dt <- fread(
  "Results/EpiMap_overlap/random.peak.freq.bed",
  col.names = c("chrom","start","end","count")
) %>%
  mutate(assay='random')

cf.dt <- fread(
  "Results/EpiMap_overlap/cfChIP.K27.peak.freq.bed",
  col.names = c("chrom","start","end","count")
) %>%
  mutate(assay='cfChIP K27')
wbc.dt <- fread(
  "Results/EpiMap_overlap/wbc.peak.freq.bed",
  col.names = c("chrom","start","end","count")
) %>%
  mutate(assay='WBC K27')

#dt = rbind.data.frame(bg.dt, cf.dt, wbc.dt)
#dt = rbind.data.frame(bg.dt, cf.dt)
dt = rbind.data.frame(wbc.dt, cf.dt)
# 2) Total number of epigenomes
N <- 98
#assay="all_downsampled"
# 3) Compute percentage detected
dt[, pct := count / N ]
dt[, length := end - start ]

# 4) Bin into categories
dt[, category := cut(
  pct,
  breaks = c(-Inf, 0, 0.25, 0.75, Inf),
  labels = c("No overlap (0%)", "Tissue-specific (1-25%)", "Moderately specific (25-75%)", "Broadly active (>75%)"),
  right = TRUE
)]

# 2) plot histogram
p_hist <- ggplot(dt%>%filter(pct>0), aes(x = pct, fill = assay)) +
  geom_histogram(
    bins       = 77,
    boundary   = 0,
    closed     = "left",
    color      = "black",
    position   = "identity",
    alpha      = 0.4
  ) +
  # add vertical lines at 25% and 75%
  geom_vline(xintercept = c(0.25, 0.75),
             color      = "red",
             linetype   = "dashed",
             size       = 1) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)  ) +
  #scale_y_log10() +
  labs(
    title = "Overlapping pct by Assay",
    x     = "Fraction of tissues active",
    y     = "Number of peaks",
    fill  = "Assay"
  ) +
  theme_classic(base_size = 14)
p_hist


dt %>%
  filter(pct!=0) %>%
  group_by(assay, category) %>%
  summarise(n(), percent())

dt %>%
  filter(pct!=0) %>%
  group_by(assay) %>%
  summarise(n())





