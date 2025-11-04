# ---- Packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(readr)
})

# ---- Working dir
setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

# ---- Inputs
peaks <- c(
  "k4_combined_cQTLs_GTEx_eqtls",
  "k27_combined_cQTLs_GTEx_eqtls",
  "wb_combined_cQTLs_GTEx_eqtls"
)

# ---- Load + tidy
dat <- lapply(peaks, function(peak) {
  fn <- file.path("Results/GTEx_eqtls_enrichment", paste0(peak, ".enrichment.summary.txt"))
  d  <- fread(fn, sep = "\t", header = FALSE,
              colClasses = c("character", rep("NULL", 8), "numeric", "NULL"))
  setnames(d, c("V1","V10"), c("Tissues","enrichment.rand"))
  d$annot <- peak
  d
}) |> rbindlist()

# tidy names
clean_names <- function(x) {
  x <- gsub("_", " ", x, fixed = TRUE)
  sub("\\sBA\\d+$", "", x)   # drop trailing " BA9" etc if present
}

dat <- dat |>
  mutate(
    tissue_cleaned = clean_names(Tissues),
    group_info = case_when(
      str_detect(tissue_cleaned, "Adipose")   ~ "Adipose",
      str_detect(tissue_cleaned, "Artery")    ~ "Artery",
      str_detect(tissue_cleaned, "Cells")     ~ "Cell line",
      str_detect(tissue_cleaned, "Brain")     ~ "Brain",
      str_detect(tissue_cleaned, "Colon")     ~ "Colon",
      str_detect(tissue_cleaned, "Esophagus") ~ "Esophagus",
      str_detect(tissue_cleaned, "Heart")     ~ "Heart",
      str_detect(tissue_cleaned, "Blood")     ~ "Whole Blood",
      str_detect(tissue_cleaned, "Skin")      ~ "Skin",
      str_detect(tissue_cleaned, "Breast")    ~ "Breast Mammary Tissue",
      TRUE ~ tissue_cleaned
    ),
    annot = factor(annot, levels = peaks,
                   labels = c("H3K4me3 cfChIP", "H3K27ac cfChIP", "H3K27ac WBC ChIP"))
  )

# mean per group × annotation; drop cell lines (per earlier spec)
dat_plot <- dat |>
  group_by(group_info, annot) |>
  summarise(value = mean(enrichment.rand, na.rm = TRUE), .groups = "drop") |>
  filter(group_info != "Cell line")

# ---- Matrix (rows = tissues, columns = annotations in fixed order)
col_order <- c("H3K27ac cfChIP", "H3K4me3 cfChIP", "H3K27ac WBC ChIP")
mat <- xtabs(value ~ group_info + annot, data = dat_plot)
mat <- mat[rev(rownames(mat)), col_order, drop = FALSE]  # reverse rows; select columns
mat_num <- as.matrix(mat)

# ---- Color scale (purple → teal → yellow), based on plotted matrix range
col_fun <- colorRamp2(
  seq(min(mat_num), max(mat_num), length.out = 100),
  colorRampPalette(c("#5e4fa2", "#66c2a5", "#fee08b"))(100)
)

# ---- In-cell labels with dynamic contrast (white on dark, black on light)
label_gp <- gpar(fontsize = 10, fontfamily = "Helvetica")
cell_fun <- function(j, i, x, y, w, h, fill) {
  # compute perceived luminance from hex "fill"
  rgb <- col2rgb(fill) / 255
  L   <- 0.2126*rgb[1] + 0.7152*rgb[2] + 0.0722*rgb[3]
  col <- ifelse(L > 0, "white", "black")
  grid.text(sprintf("%.1f", mat_num[i, j]), x, y, gp = gpar(col = col, fontsize = 10, fontfamily = "Helvetica"))
}

ht <- Heatmap(
  mat_num,
  name = "Fold enrichment",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 0,
  column_names_centered = TRUE,
  cell_fun = cell_fun,
  border = TRUE,
  # <- use these:
  row_names_gp    = gpar(fontsize = 10, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
  row_title_gp    = gpar(fontsize = 10, fontfamily = "Helvetica"),
  column_title_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
  heatmap_legend_param = list(
    title = "Fold enrichment \nat eQTLs",
    title_gp  = gpar(fontsize = 10, fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 9,  fontfamily = "Helvetica")
  )
)

ht
# ---- Save vector PDF (journal-friendly). Size ~120×120 mm.
pdf("Manuscripts/Figures/3_e.pdf", width = 180/25.4, height = 120/25.4, family = "Helvetica", useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right")
dev.off()

# ---- Export table (wide) for supplement
dat_wide <- dat_plot |>
  tidyr::pivot_wider(names_from = annot, values_from = value)
write_excel_csv(dat_wide, "Tables/sup7.csv")



library(ComplexHeatmap)
library(grid)

col_avg <- dat %>%
  group_by(annot) %>%
  summarise(value = round(mean(enrichment.rand, na.rm = TRUE), digits = 2))
col_avg = col_avg$value[c(2,1,3)]
top_anno <- HeatmapAnnotation(
  `Mean ` = anno_barplot(
    col_avg,
    ylim = c(0, max(col_avg, na.rm = TRUE)),
    gp = gpar(fill = "skyblue4", col = NA),
    bar_width = 0.6
  ),
  # optional: also print the numbers above each bar
  Mean = anno_text(sprintf("%.2f", col_avg), rot = 0,
                   gp = gpar(fontsize = 9, fontfamily = "Helvetica")),
  annotation_name_side = "left",
  annotation_height = unit(c(10, 6), "mm")
)

ht <- Heatmap(
  mat_num,
  name = "Fold enrichment",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = top_anno,
  row_names_side = "left",
  column_names_rot = 0,
  column_names_centered = TRUE,
  cell_fun = cell_fun,
  border = TRUE,
  # <- use these:
  row_names_gp    = gpar(fontsize = 10, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
  row_title_gp    = gpar(fontsize = 10, fontfamily = "Helvetica"),
  column_title_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
  heatmap_legend_param = list(
    title = "Fold enrichment \nat eQTLs",
    title_gp  = gpar(fontsize = 10, fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 9,  fontfamily = "Helvetica")
  )
)

ht
pdf("Manuscripts/Figures/3_e_with_MeanBar.pdf", width = 180/25.4, height = 140/25.4, family = "Helvetica", useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right")
dev.off()
  
