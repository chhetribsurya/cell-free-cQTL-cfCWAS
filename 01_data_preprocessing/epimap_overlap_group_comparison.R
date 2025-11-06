setwd("~/Projects/cfChIP")
# visualize_detection.R
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)   # for percent()
library(RColorBrewer)
library(tidyverse)
library(nnet)
library(ggpubr)
library(extrafont)  # For managing fonts
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 9)
)
group_map = fread('Results/EpiMap_overlap/back_up/EpiMap_info.csv')[,1:4]
colnames(group_map)=c("EID", "NAME", "GROUP", "Blood")
N=98
threshold=0.25
# 1) Read in the counts
wb.dt <- fread(
  "Results/EpiMap_overlap/wbc.peak.overlaped.tsv"
)%>%
  mutate(row_sum=rowSums(across(4:101)),
         pct=row_sum/N,
         assay='WBC K27')

cf.dt <- fread(
  "Results/EpiMap_overlap/cfChIP_H3K27ac.peak.overlaped.tsv"
) %>%
  mutate(row_sum=rowSums(across(4:101)),
         pct=row_sum/N,
         assay='cfChIP K27')

wb.top <- wb.dt %>%
  filter(pct > 0) %>%
  arrange(desc(pct)) %>%
  slice_tail(prop = threshold)

cf.top <- cf.dt %>%
  filter(pct > 0) %>%
  arrange(desc(pct)) %>%
  slice_tail(prop = threshold)

# Count number of peaks detected per EID in top 10%
counts.wb <- apply(wb.top[, 4:101], 2, sum)
counts.cf <- apply(cf.top[, 4:101], 2, sum)

eid <- colnames(wb.dt)[4:101]

df <- data.frame(
  EID = eid,
  wb = as.numeric(counts.wb),
  cf = as.numeric(counts.cf)
)

df <- left_join(df, group_map, by = "EID")
df_long <- df %>%
  pivot_longer(cols = c("wb", "cf"), names_to = "source", values_to = "overlap")

blood.table = df_long %>%
  group_by(Blood, source) %>%
  summarise(total = mean(overlap))


df.group = df %>%
  group_by(Blood) %>%
  summarise(
    mean_diff = mean(cf - wb),
    p_value = if (n() > 1) wilcox.test(cf, wb)$p.value else NA_real_
  )%>%
  mutate(label = paste0("Mean diff. = ", round(mean_diff, 1), 
                        "\nP = ", format.pval(p_value, digits = 2, eps = 1e-30)))

p <- ggplot(aes(x = Blood, y = total, fill = source), data=blood.table) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    name   = "Assays",
    values = c("cf" = "#1f78b4",
               "wb" = "#e31a1c"),
    labels = c("cf" = "cfChIP", "wb" = "WBC ChIP")
  )+
  scale_x_discrete(labels = c("B" = "Blood cells", "NB" = "Non-blood tissues"))+
  labs(
    x  = "Category",
    y = "The mean number of chromatin peaks overlapped"
  )+
  theme(
    # We keep font sizes at 7 (or smaller if needed)
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.text.y  = element_text(size = 9)
  )+
  geom_text(
    data = df.group,
    aes(x = Blood, y = Inf, label = label),
    vjust = 1,    # adjust label height above the tallest bar
    hjust = 0.3,
    size = 3,
    inherit.aes = FALSE
  )
p

ggsave(
  filename = paste0('Figures/EpiMap Overlap/NB_B_group_overlapped_mean_', threshold, '.pdf'),
  plot = p,
  device = pdf,  # Vector output
  units = "mm",
  width = 110,
  height = 160
)



counts.wb <- apply(wb.dt[, 4:101], 2, sum)
counts.cf <- apply(cf.dt[, 4:101], 2, sum)

eid <- colnames(wb.dt)[4:101]

df <- data.frame(
  EID = eid,
  wb = as.numeric(counts.wb),
  cf = as.numeric(counts.cf)
)

df <- left_join(df, group_map, by = "EID")
