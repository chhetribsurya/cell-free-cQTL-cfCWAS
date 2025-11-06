library(ggplot2)
library(readr)   # if you print ratios etc.
setwd("~/Projects/cfChIP")
# Global text settings
base_family <- "Helvetica"
base_size   <- 10
fig_w_mm    <- 75
fig_h_mm    <- 75

theme_set(theme_void(base_family = base_family, base_size = base_size))

make_pie <- function(df, fills, file_pdf, show_legend = FALSE, title = NULL) {
  p <- ggplot(df, aes(x = "", y = number, fill = label)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = fills, name = NULL) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = base_size, colour = "black"),
      plot.title  = element_text(size = base_size, colour = "black")
    ) +
    labs(title = title)
  
  ggsave(
    filename = file_pdf,
    plot = p,
    device = "pdf",            # native PDF on macOS; no XQuartz required
    family = base_family,
    useDingbats = FALSE,
    units = "mm",
    width = fig_w_mm, height = fig_h_mm
  )
  p
}

## ---- Pie 1: H3K27ac peaks ----
df <- data.frame(
  Peaks  = c("H3K27ac peaks", "Imbalanced/cQTL peak"),
  number = c(51542, 4816)
)
df$label <- paste0(df$Peaks, ": ", df$number)
print(df$number[2] / df$number[1])

make_pie(
  df,
  fills = c("H3K27ac peaks: 51542" = "#999999",
            "Imbalanced/cQTL peak: 4816" = "#e9c716"),
  file_pdf = "Manuscripts/Figures/2_b_pie_H3K27ac.pdf",
  show_legend = FALSE
)

## ---- Pie 2: H3K4me3 peaks ----
df <- data.frame(
  Peaks  = c("H3K4me3 peaks", "Imbalanced/cQTL peak"),
  number = c(51888, 7357)
)
print(df$number[2] / df$number[1])
df$label <- paste0(df$Peaks, ": ", df$number)

make_pie(
  df,
  fills = c("H3K4me3 peaks: 51888" = "#999999",
            "Imbalanced/cQTL peak: 7357" = "#50ad9f"),
  file_pdf = "Manuscripts/Figures/2_b_pie_H3K4me3.pdf",
  show_legend = FALSE
)

