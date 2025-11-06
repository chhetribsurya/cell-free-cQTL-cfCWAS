# --- Setup
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(purrr)
})

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))
mm <- function(x) x/25.4

setwd("~/Projects/cfChIP")
cwas_dir <- "Results/CWAS/cwas_resuts_peak/"

blood_cell.list <- c("30000_irnt","30120_irnt",
                 "30130_irnt","30140_irnt")


panel.list <- list.files(cwas_dir)

# --- Small safe divide helper
percent <- function(sig, total) if (total > 0) (sig/total)*100 else 0

# --- Loader → tall summary (panel, trait, no.cv, no.sig, prop.sig, prop.sig.permile)
summarise_set <- function(trait_vec) {
  map_dfr(panel.list, function(panel) {
    map_dfr(trait_vec, function(trait) {
      sig.fp <- file.path(cwas_dir, panel, paste0(trait, ".cwas.sig.best.txt"))
      cv.fp  <- file.path(cwas_dir, panel, paste0(trait, ".cwas.cv.best.txt"))
      if (!file.exists(sig.fp) || !file.exists(cv.fp)) return(NULL)
      sig.df <- read.table(sig.fp, header = TRUE)
      cv.df  <- read.table(cv.fp,  header = TRUE)
      no.sig <- nrow(sig.df); no.cv <- nrow(cv.df)
      tibble(
        panel = panel,
        trait = trait,
        no.cv = no.cv,
        no.sig = no.sig,
        prop.sig = if (no.cv > 0) no.sig/no.cv else 0,
        prop.sig.percent = percent(no.sig, no.cv)
      )
    })
  })
}

blood_cell.df    <- summarise_set(blood_cell.list)

# --- Nice labels
pretty_trait <- function(x) dplyr::case_when(
  x == "30000_irnt" ~ "Leukocyte count",
  x == "30120_irnt" ~ "Lymphocyte count",
  x == "30130_irnt" ~ "Monocyte count",
  x == "30140_irnt" ~ "Neutrophill count",
  TRUE ~ x
)

pretty_panel <- function(x) dplyr::case_when(
  x == "cfChIP_H3K27ac" ~ "cfChIP",
  x == "whole-blood"    ~ "WBC ChIP",
  TRUE ~ x
)

# --- Plot helper
make_bar <- function(df, outfile_pdf) {
  plot.df <- df %>%
    filter(panel %in% c("cfChIP_H3K27ac","whole-blood")) %>%
    mutate(
      panel = pretty_panel(panel),
      trait = pretty_trait(trait)
    )
  
  # order traits by cfChIP permille (ascending, then flip for descending if you want)
  order_df <- plot.df %>%
    filter(panel == "cfChIP") %>%
    arrange(prop.sig.percent)
  trait_levels <- rev(order_df$trait)  # highest at top/right
  plot.df <- plot.df %>% mutate(trait = factor(trait, levels = trait_levels))
  
  p <- ggplot(plot.df, aes(x = trait, y = prop.sig.percent, fill = panel)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.8, colour = "black", linewidth = 0.2) +
    scale_fill_manual(
      values = c("cfChIP" = "#e9c716", "WBC ChIP" = "firebrick"),
      breaks = c("cfChIP","WBC ChIP"),
      name = ""
    ) +
    labs(
      x = "GWAS",
      y = "Significant CWAS associations \nper 100 genetic models"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title  = element_text(size = 10, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 10, colour = "black"),
      legend.text  = element_text(size = 10, colour = "black"),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  ggsave(outfile_pdf, p, device = "pdf", width = mm(150), height = mm(120), dpi = 300)
  print(p)
}

# --- Build and save
p_blood    <- make_bar(blood_cell.df,    "Manuscripts/Figures/6_a_blood.pdf")
