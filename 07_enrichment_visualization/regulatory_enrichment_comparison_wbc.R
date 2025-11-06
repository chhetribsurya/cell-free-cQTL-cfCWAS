# Set working directory
setwd("~/Projects/cfChIP")

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
  library(purrr)
  library(extrafont)
})
source('../Scripts/functions.R')

# Use a theme with Helvetica, base size of 7pt
theme_set(theme_classic(base_family = "Helvetica", base_size = 7))

# --- Inputs
cancer.list <- c(
  "PASS_OvarianCancer",
  "ProstateCancer_Meta_Schumacher2018.nodup",
  "UKB_460K.cancer_MELANOMA",
  "breast_cancer",
  "ec",
  "lc",
  "rcc_meta"
)

noncancer.list <- c("PASS_Alzheimers_Jansen2019","PASS_Crohns_Disease","PASS_IBD_deLange2017",
                    "PASS_Rheumatoid_Arthritis","PASS_Schizophrenia","PASS_Ulcerative_Colitis",
                    "UKB_460K.disease_PSORIASIS")

cwas_dir <- "Results/CWAS/cwas_resuts_peak/"
gwas_dir <- "Results/CWAS/GWAS_sig_loci/"

# --- Read CWAS cancer results (all panels)
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

# --- Read GWAS cancer results
cancer.gwas.df <- map_dfr(cancer.list, function(cancer) {
  fp <- file.path(gwas_dir, paste0(cancer, ".sig.snps.bed"))
  if (file.exists(fp)) {
    read.table(fp, header = FALSE, col.names = c("CHR","P0","P1")) %>%
      mutate(trait = cancer, panel = "GWAS risk loci")
  } else NULL
})

# --- Combine minimal columns (no cytobands)
df <- bind_rows(
  cancer.df %>% select(CHR, P0, P1, trait, panel),
  cancer.gwas.df %>% select(CHR, P0, P1, trait, panel)
)

# --- Pretty trait names
rename.list <- c("Ovarian Cancer","Prostate Cancer","Melanoma","Breast Cancer",
                 "Endometrial Cancer","Lung Cancer","Renal Cell Carcinoma")
cancer_rename_map <- setNames(rename.list, cancer.list)
df <- df %>% mutate(trait = recode(trait, !!!cancer_rename_map))

# --- Helper: merge intervals within ±500 kb and keep original peak positions
merge_within_pm500kb <- function(df_in, flank_bp = 500000L) {
  if (nrow(df_in) == 0) {
    return(tibble(trait = character(), panel = character(),
                  locus_id = character(), CHR = character(), P0 = integer(), P1 = integer()))
  }
  # Build GRanges with explicit columns
  gr <- makeGRangesFromDataFrame(df_in,
                                 seqnames.field = "CHR",
                                 start.field    = "P0",
                                 end.field      = "P1",
                                 keep.extra.columns = TRUE,
                                 ignore.strand  = TRUE)
  # Expand each interval by ±flank
  start(gr) <- pmax(1L, start(gr) - flank_bp)
  end(gr)   <- end(gr) + flank_bp
  
  # Reduce with revmap to know which original peaks belong to each merged locus
  grm <- GenomicRanges::reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
  locus_id <- paste0(as.character(seqnames(grm)), ":", start(grm), "-", end(grm))
  
  # Map back to original rows; keep the original CHR/P0/P1
  revmap_list <- S4Vectors::mcols(grm)$revmap
  out <- map2_dfr(revmap_list, locus_id, function(ix, id) {
    orig <- as.data.frame(S4Vectors::mcols(gr)[ix, c("trait","panel")])
    # original positions from the unexpanded intervals:
    orig_pos <- as.data.frame(gr)[ix, c("seqnames","start","end")]
    tibble(
      trait   = as.character(orig$trait),
      panel   = as.character(orig$panel),
      locus_id = id,
      CHR     = as.character(orig_pos$seqnames),
      P0      = as.integer(orig_pos$start),
      P1      = as.integer(orig_pos$end)
    )
  })
  out
}

# --- Remove cfMeDIP for downstream plots
df_no_medip <- df %>% filter(panel != "cfMeDIP")

# Choose the cancer to visualize in the Venn
cancer_to_plot <- "Prostate Cancer"

# Subset to that cancer
sub_df <- df_no_medip %>% filter(trait == cancer_to_plot)

# Compute merged loci (±500 kb) and retain original peak info
merged_map <- merge_within_pm500kb(sub_df, flank_bp = 500000L)
# merged_map has: trait, panel, locus_id, CHR, P0, P1 (original peak positions)

# Build sets for Venn by merged locus_id (three panels: cfChIP H3K27ac, cfChIP H3K4me3, GWAS)
panel_name_map <- c(
  "cfChIP_H3K27ac" = "cfChIP H3K27ac",
  "cfChIP_H3K4me3" = "cfChIP H3K4me3",
  "GWAS risk loci" = "GWAS"
)
merged_map <- merged_map %>% mutate(panel = recode(panel, !!!panel_name_map))

venn_sets <- list(
  `cfChIP H3K27ac` = unique(merged_map$locus_id[merged_map$panel == "cfChIP H3K27ac"]),
  `cfChIP H3K4me3` = unique(merged_map$locus_id[merged_map$panel == "cfChIP H3K4me3"]),
  GWAS             = unique(merged_map$locus_id[merged_map$panel == "GWAS"])
)
# Get uniquely present elements for each set
unique_elems <- lapply(names(venn_sets), function(set_name) {
  others <- unlist(venn_sets[names(venn_sets) != set_name])
  setdiff(venn_sets[[set_name]], others)
})

names(unique_elems) <- names(venn_sets)
unique_elems
# Draw Venn (3 sets; scientific-friendly colors)
venn_plot <- ggvenn(
  venn_sets,
  fill_color   = c("#e9c716", "#50ad9f", "#377EB8"),
  stroke_color = "black",
  stroke_size  = 0.5,
  set_name_size = 3.5,
  text_size     = 3.5,
  show_percentage = FALSE
)

print(venn_plot)

# Save as vector PDF (120 x 100 mm if you prefer that size; you had 150 x 120 mm before)
ggsave("Manuscripts/Figures/4_e.pdf",
       plot = venn_plot,
       device = pdf,
       width = 100, height = 100, units = "mm", dpi = 300)

# ----- OPTIONAL: export the mapping table for supplement
# write.csv(merged_map, "Tables/sup_locus_mapping_pm500kb.csv", row.names = FALSE)
