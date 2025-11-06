# Set working directory
setwd("~/Projects/cfChIP")

# Libraries
library(biovizBase)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(purrr)
source('../Scripts/functions.R')
library(extrafont)

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

# Traits
cancer.list <- c(
  "PASS_OvarianCancer", "ProstateCancer_Meta_Schumacher2018.nodup",
  "UKB_460K.cancer_MELANOMA", "breast_cancer", "ec", "lc", "rcc_meta"
)

# Directories
cwas_dir <- "Results/CWAS/cwas_resuts_peak/"
gwas_dir <- "Results/CWAS/GWAS_sig_loci/"

# -------------------------------
# Function to merge within window
# -------------------------------
merge_within_window <- function(df, window_size = 1e6) {
  # Rename for GRanges
  
  gr <- makeGRangesFromDataFrame(
    df,
    seqnames.field    = "CHR",
    start.field       = "P0",
    end.field         = "P1",
    keep.extra.columns = TRUE,
    ignore.strand     = TRUE
  )
  gr_merged <- GenomicRanges::reduce(gr, min.gapwidth = window_size, with.revmap = TRUE)
  
  merged_df <- as.data.frame(gr_merged) %>%
    mutate(locus_id = paste0(seqnames, ":", start, "-", end)) %>%
    select(locus_id)
  
  revmap_list <- mcols(gr_merged)$revmap
  trait_panel <- lapply(seq_along(revmap_list), function(i) {
    orig_rows <- df[revmap_list[[i]], c("trait", "panel")]
    orig_rows$locus_id <- merged_df$locus_id[i]
    return(orig_rows)
  }) %>% bind_rows()
  
  return(trait_panel)
}

# -------------------------------
# Load CWAS cancer results
# -------------------------------
panel.list <- list.files(cwas_dir)
cancer.df <- map_dfr(panel.list, function(panel) {
  map_dfr(cancer.list, function(cancer) {
    fp <- file.path(cwas_dir, panel, paste0(cancer, ".cwas.sig.best.txt"))
    if (file.exists(fp)) {
      read.table(fp, header = TRUE) %>%
        mutate(panel = panel, trait = cancer)
    } else NULL
  })
})

# Load GWAS cancer results
cancer.gwas.df <- map_dfr(cancer.list, function(cancer) {
  fp <- file.path(gwas_dir, paste0(cancer, ".sig.snps.bed"))
  if (file.exists(fp)) {
    read.table(fp, header = FALSE, col.names = c("CHR","P0","P1")) %>%
      mutate(trait = cancer)
  } else NULL
})

# -------------------------------
# Merge to 1 Mb loci
# -------------------------------
cancer.gwas.df <- cancer.gwas.df %>% mutate(panel = "GWAS risk loci")
gwas_loci <- merge_within_window(cancer.gwas.df, window_size = 1e6)

# Remove cfMeDIP here
cancer.df <- cancer.df %>% filter(panel != "cfMeDIP", panel != "whole-blood")
cwas_loci <- merge_within_window(cancer.df, window_size = 1e6)

# -------------------------------
# Novelty classification
# -------------------------------
novel.summary <- cancer.df %>%
  group_by(panel, trait, CHR, P0, P1) %>%
  summarise(max_gwas_abs_z = max(abs(BEST.GWAS.Z), na.rm = TRUE), .groups = "drop") %>%
  mutate(novel = if_else(max_gwas_abs_z < 5.45, "GWAS sub-significant", "GWAS significant"))

# -------------------------------
# Pretty names
# -------------------------------
panel_map <- c("cfChIP_H3K27ac"="cfChIP H3K27ac",
               "cfChIP_H3K4me3"="cfChIP H3K4me3")
trait_map <- c(
  "PASS_OvarianCancer"="Ovarian Cancer",
  "ProstateCancer_Meta_Schumacher2018.nodup"="Prostate Cancer",
  "UKB_460K.cancer_MELANOMA"="Melanoma",
  "breast_cancer"="Breast Cancer",
  "ec"="Endometrial Cancer",
  "lc"="Lung Cancer",
  "rcc_meta"="Renal Cell Carcinoma"
)

novel.summary <- novel.summary %>%
  mutate(panel = recode(panel, !!!panel_map),
         trait = recode(trait, !!!trait_map))

# -------------------------------
# Counts per trait/panel/novelty
# -------------------------------
novel.plot <- novel.summary %>%
  group_by(panel, trait, novel) %>%
  summarise(N = n(), .groups = "drop")

# Order traits by cfChIP H3K27ac total
order_vec <- novel.plot %>%
  filter(panel == "cfChIP H3K27ac") %>%
  group_by(trait) %>%
  summarise(total = sum(N), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(trait)

novel.plot <- novel.plot %>%
  mutate(trait = factor(trait, levels = order_vec))

# -------------------------------
# Helper plotting function
# -------------------------------
plot_panel <- function(df_panel, sig_col = "#377EB8") {
  ggplot(df_panel, aes(x = trait, y = N, fill = novel)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("GWAS sub-significant" = sig_col,
                                 "GWAS significant"     = "grey70"),
                      breaks = c("GWAS significant","GWAS sub-significant"),
                      name = NULL) +
    labs(x = NULL, y = "Loci count") +
    theme_classic(base_family = "Helvetica", base_size = 10) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title  = element_text(size = 10, colour = "black"),
      legend.text = element_text(size = 10, colour = "black"),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm")
    )
}

# -------------------------------
# Plots & save (size in mm)
# -------------------------------
mm_to_in <- function(mm) mm / 25.4

# cfChIP H3K27ac
p1 <- novel.plot %>% filter(panel == "cfChIP H3K27ac") %>% plot_panel()
p1
ggsave("Manuscripts/Figures/ext_3_c.pdf", p1, device = pdf,
       width = mm_to_in(120), height = mm_to_in(80),
       family = "Helvetica", useDingbats = FALSE,)

# cfChIP H3K4me3
p2 <- novel.plot %>% filter(panel == "cfChIP H3K4me3") %>% plot_panel()
p2
ggsave("Manuscripts/Figures/4_d.pdf", p2, device = pdf,
       width = mm_to_in(120), height = mm_to_in(100), dpi = 300)
