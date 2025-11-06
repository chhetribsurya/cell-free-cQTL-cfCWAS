# Set working directory
setwd("~/Projects/cfChIP")

# Load required libraries
library(biovizBase)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(purrr)  
source('../Scripts/functions.R')
library(extrafont)  

# Set theme: Helvetica, base size 10pt
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 10)
)

# Define lists of cancer traits
cancer.list <- c(
  "PASS_OvarianCancer",
  "ProstateCancer_Meta_Schumacher2018.nodup",
  "UKB_460K.cancer_MELANOMA",
  "breast_cancer",
  "ec",
  "lc",
  "rcc_meta"
)

# Define base directories
cwas_dir <- "Results/CWAS/cwas_resuts_peak/"
gwas_dir <- "Results/CWAS/GWAS_sig_loci/"

# Read CWAS cancer results
panel.list <- list.files(cwas_dir)
cancer.df <- map_dfr(panel.list, function(panel) {
  map_dfr(cancer.list, function(cancer) {
    file_path <- file.path(cwas_dir, panel, paste0(cancer, ".cwas.sig.best.txt"))
    if (file.exists(file_path)) {
      df <- read.table(file_path, header = TRUE)
      df <- df %>% mutate(panel = panel, trait = cancer)
      return(df)
    } else {
      return(NULL)
    }
  })
})

# Read GWAS cancer results
cancer.gwas.df <- map_dfr(cancer.list, function(cancer) {
  file_path <- file.path(gwas_dir, paste0(cancer, ".sig.snps.bed"))
  if (file.exists(file_path)) {
    df <- read.table(file_path, header = FALSE, col.names = c("CHR", "P0", "P1"))
    df <- df %>% mutate(trait = cancer)
    return(df)
  } else {
    return(NULL)
  }
})

# Merge loci within ±500 kb (1 Mb window total)
# Merge signals into loci if they are within ±500 kb (i.e., <= 1 Mb apart)
merge_within_window <- function(df, flank_bp = 500000L) {
  if (is.null(df) || nrow(df) == 0) {
    return(dplyr::tibble(trait = character(), panel = character(), locus = character()))
  }
  # Build GRanges explicitly naming the columns
  gr <- makeGRangesFromDataFrame(
    df,
    seqnames.field    = "CHR",
    start.field       = "P0",
    end.field         = "P1",
    keep.extra.columns = TRUE,
    ignore.strand     = TRUE
  )
  
  # Expand each interval by ±flank_bp
  start(gr) <- pmax(1L, start(gr) - flank_bp)
  end(gr)   <- end(gr) + flank_bp
  
  # Merge overlapping/nearby expanded intervals; keep mapping back to originals
  grm <- GenomicRanges::reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
  
  # Locus IDs like "chr1:12345-67890"
  locus_id <- paste0(as.character(seqnames(grm)), ":", start(grm), "-", end(grm))
  
  # Map merged loci back to their original trait/panel rows
  revmap_list <- S4Vectors::mcols(grm)$revmap
  out <- purrr::map2_dfr(revmap_list, locus_id, function(ix, id) {
    dplyr::tibble(
      trait = as.character(S4Vectors::mcols(gr)$trait[ix]),
      panel = as.character(S4Vectors::mcols(gr)$panel[ix]),
      locus = id
    )
  })
  
  out
}


# Process GWAS loci
cancer.gwas.df <- cancer.gwas.df %>% mutate(panel = "GWAS risk loci")
cancer.gwas.df <- merge_within_window(cancer.gwas.df)

# Process CWAS loci
cancer.df <- merge_within_window(cancer.df)

# Combine all panels
df <- bind_rows(cancer.df, cancer.gwas.df)

# Rename traits
rename.list <- c(
  "Ovarian Cancer",
  "Prostate Cancer",
  "Melanoma",
  "Breast Cancer",
  "Endometrial Cancer",
  "Lung Cancer",
  "Renal Cell Carcinoma"
)
cancer_rename_map <- setNames(rename.list, cancer.list)
df <- df %>% mutate(trait = recode(trait, !!!cancer_rename_map))

# Summarize unique loci counts
df_summary <- df %>%
  group_by(panel, trait) %>%
  summarise(unique_cytoband_count = n_distinct(locus), .groups = "drop")

# Order traits by GWAS loci counts
gwas_cytoband_counts <- df_summary %>%
  filter(panel == "GWAS risk loci") %>%
  arrange(desc(unique_cytoband_count)) %>%
  pull(trait)
df_summary$trait <- factor(df_summary$trait, levels = gwas_cytoband_counts)

# Define panel order and remove cfMeDIP + WBC H3K27ac
panel_levels <- c("GWAS risk loci", "cfChIP_H3K4me3", "cfChIP_H3K27ac")
panel_rename_map <- c(
  "GWAS risk loci" = "GWAS risk loci",
  "cfChIP_H3K4me3" = "cfChIP H3K4me3",
  "cfChIP_H3K27ac" = "cfChIP H3K27ac"
)
df_summary$panel <- factor(df_summary$panel, levels = panel_levels)
df_summary = df_summary[!is.na(df_summary$panel),]
df_summary$panel <- recode(df_summary$panel, !!!panel_rename_map)

# Plot
p <- ggplot(df_summary, aes(x = trait, y = unique_cytoband_count, fill = panel)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  scale_fill_manual(values = c("#377EB8", "#50ad9f", "#e9c716")) +
  labs(y = "Loci counts (±500 kb)", fill = "Panel") +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.ticks = element_line(),
    legend.key.size = unit(0.3, "cm")
  )
p
# Save
ggsave("Manuscripts/Figures/4_c.pdf", plot = p,
       device = "pdf", family = "Helvetica", useDingbats = FALSE,
       units  = "mm", width = 120, height = 80)   # adjust size as needed)
