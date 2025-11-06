library(tidyverse)
library(data.table)

# ---- Theme: Helvetica, 10 pt (journal-friendly)
theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

# ---- Set working directory
setwd("~/Projects/cfChIP")

# ---- Input
# File contains: Roadmap group name, % overlap (Measurement), and panel ("cfChIP_H3K27ac" or "whole-blood")
res.file <- "Results/EpiMap_overlap/EnhA1_EnhA2_EnhG1_EnhG2_grouped_output_ratios.txt"

# ---- Load data
res <- fread(res.file, col.names = c("Group", "Measurement", "Panel")) %>%
  as_tibble()

res  = res %>%
  filter(Group != "ENCODE cell lines")
# ---- Reshape to wide format
df_wide <- res %>%
  pivot_wider(names_from = Panel, values_from = Measurement) %>%
  mutate(
    ratio     = cfChIP_H3K27ac / `whole-blood`,
    logratio  = log(ratio),
    diff      = (cfChIP_H3K27ac - `whole-blood`)*100,
    sign      = factor(if_else(diff > 0, "cfChIP > WBC", "cfChIP < WBC"),
                       levels = c("cfChIP > WBC", "cfChIP < WBC"))
  ) %>%
  arrange(diff) %>%
  mutate(Group = factor(Group, levels = Group)) # order by diff

# ---- Plot
p <- ggplot(df_wide, aes(x = Group, y = diff, fill = sign)) +
  geom_col(width = 0.75, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "black") +
  coord_flip() +
  scale_fill_manual(
    values = c("cfChIP > WBC" = "#1f78b4",
               "cfChIP < WBC" = "#e31a1c"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(-8, 11), expand = expansion(mult = c(0, 0))) +
  labs(
    y = expression(
      atop(
        "Difference in % cQTL peaks overlapped with merged chromatin",
        "(" * cfChIP - WBC * ")"
      )
    ),
    x = "EpiMap groups"
  )+
  theme(
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x  = element_text(size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    legend.position = "bottom",
    legend.text  = element_text(size = 10, colour = "black"),
    plot.margin  = margin(5, 10, 5, 5)
  )

p

# ---- Save vector PDF
ggsave(
  filename = "Manuscripts/Figures/3_a.pdf",
  plot = p,
  device = "pdf",
  family = "Helvetica",
  useDingbats = FALSE,
  units = "mm",
  width = 120, height = 120
)
