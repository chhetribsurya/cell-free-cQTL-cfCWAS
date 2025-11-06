library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(readxl)
library(ggrepel)
library(extrafont)  # To manage fonts if needed

# Use a theme with Helvetica, base size of 7pt
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 9)
)

setwd("~/Projects/cfChIP")

state <- 'EnhA1_EnhA2_TssA_EnhG1_EnhG2'
result.file <- paste0('~/Projects/cfChIP/Results/EpiMap_overlap/',
                      state, '_combined_results.csv')
df <- read.csv(result.file, header = TRUE)

df <- df %>%
  mutate(CF_freq = CF_count / Total_count,
         WB_freq = WB_count / Total_count,
         ratio   = CF_count / WB_count,
         delta   = (CF_freq - WB_freq) * 100,
         logratio = log(CF_count / WB_count))

##############
# Create the bar plot
pie_df <- data.frame(
  category = c("cfChIP > WBC ChIP", "WBC ChIP > cfChIP"),
  lab = c("cfChIP-enriched tissues", "WBC-enriched tissues "),
  count    = c(sum(df$delta>0), sum(df$delta<0))
) %>%
  mutate(pct   = count / sum(count),
         label = sprintf("%s (%s)", lab, scales::percent(pct, 1)),
         mid   = cumsum(pct) - pct/2)   # slice midpoints for polar y

cut <- 0.08  # threshold for "small" slices (8%)

p <- ggplot(pie_df, aes(x = 1, y = pct, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +
  
  # ----- labels for larger slices (centered INSIDE) -----
geom_segment(
  data = subset(pie_df, pct >= cut),
  aes(x = 1.7, xend = 1.3, y = 0.5, yend = 0.5),
  inherit.aes = FALSE,
  color="black"
) +
geom_label(
  data = subset(pie_df, pct >= cut),
  aes(x=1.7, y = 0.8, label = label),
  position = position_stack(vjust = 0.6),
  size = 10/.pt,
  label.size = NA,                    # border width
  label.r = grid::unit(0, "pt"),       # rounded corner radius
  label.padding = grid::unit(4, "pt"), # same padding for all boxes
  fill = "white", color = "black"
) +
  
  # ----- labels for tiny slices (OUTSIDE at fixed radius) -----
# leader line
geom_segment(
  data = subset(pie_df, pct < cut),
  aes(x = 1.7, xend = 1.3, y = 0.025, yend = 0.025),
  inherit.aes = FALSE,
  color="black"
) +
  # boxed label at x = 1.28 (outside the pie)
  geom_label(
    data = subset(pie_df, pct < cut),
    aes(x = 1.7, y = mid, label = label),
    size = 10/.pt,
    label.size = NA,
    label.r = grid::unit(0, "pt"),
    label.padding = grid::unit(4, "pt"),
    fill = "white", color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = c("cfChIP > WBC ChIP" = "#1f78b4",
                               "WBC ChIP > cfChIP" = "#e31a1c"),
                    name = NULL) +
  expand_limits(x = 1.45) +   # room for the outside label
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size=10))

p

ggsave(
  filename = paste0('Manuscripts/Figures/1_d.pdf'),
  plot = p,
  device = pdf,  # Vector output
  units = "mm",
  family = "Helvetica",      # use mac Helvetica
  useDingbats = FALSE,       # avoids weird font substitutions
  width = 100,
  height = 100
)

