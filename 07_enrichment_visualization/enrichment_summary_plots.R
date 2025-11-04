# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

regulatory_element <- "active"
regulatory_element <- "enhancer"
regulatory_element <- "promoter"
regulatory_element <- "bivalent"

base_dir <- "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cQTL_H3K4me3_qvalBased_V2"
base_dir <- "/Users/chhetribsurya/Dropbox/github_repo/dfci-harvard/cell-free_cQTL-CWAS/NewTest_cfChIP-WBC_peaks_H3K4me3Enrich_V2"
# base_dir <- "/project/NewTest_cfChIP-WBC_peaks_H3K4me3Enrich"

# Read the data (only once, use the full random stats file)
data <- read.csv(file.path(base_dir, regulatory_element, "detailed_statistics_full_random.csv"))

# Clean category names
data$category <- gsub("_", " ", data$category)

# Fixed color mapping for each category
category_colors <- c(
  "fetal specific" = "#F1C40F",
  "developmental specific" = "#E74C3C",
  "stem cell specific" = "#3498DB",
  "adult specific" = "#2ECC71"
)

# Ensure the factor levels match the color mapping order
data$category <- factor(data$category, levels = names(category_colors))

# Set output directory
output_dir <- file.path(base_dir, regulatory_element, "plots", "enrich")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Fold Enrichment Plot
# Sort data by fold_enrichment
data_fold <- data[order(data$fold_enrichment),]

pdf(file.path(output_dir, "fold_enrichment_plot.pdf"), width=16, height=14)
    par(mar=c(10,4,4,6))  # Increased bottom margin for x-axis labels
    n <- nrow(data_fold)
    bp <- barplot(data_fold$fold_enrichment, 
                  names.arg=rep('', n), # suppress default labels
                  col=category_colors[as.character(data_fold$category)],
                  width=0.05,   # slim bars
                  space=0.2,   # reduced space between bars
                  ylim=c(0, max(data_fold$fold_enrichment) * 1.2),
                  main="Fold Enrichment by Category",
                  xlab="",
                  ylab="Fold Enrichment",
                  axes=FALSE)
    axis(2)
    abline(h=0, lwd=1)
    abline(h=1, lty=2, col="red", lwd=2)

    # Add x-axis ticks at bar centers
    axis(1, at=bp, labels=FALSE, tck=-0.02)

    # Add 45-degree rotated x-axis labels at bar centers
    text(x=bp, y=par("usr")[3] - 0.05 * max(data_fold$fold_enrichment),
         labels=data_fold$category, srt=45, adj=1, xpd=TRUE, cex=1.2)

    # Add value labels
    text(x=bp, y=data_fold$fold_enrichment + 0.1,
         labels=round(data_fold$fold_enrichment, 2),
         pos=3, cex=1.5)

    # Add legend
    legend("topleft",
           legend=names(category_colors),
           fill=category_colors,
           cex=1.5,
           bty="n",
           xpd=TRUE)
dev.off()

# 2. Odds Ratio Plot
# Sort data by odds_ratio
data_odds <- data[order(data$odds_ratio),]

pdf(file.path(output_dir, "odds_ratio_plot.pdf"), width=16, height=14)
    par(mar=c(10,4,4,6))
    n <- nrow(data_odds)
    bp <- barplot(data_odds$odds_ratio, 
                  names.arg=rep('', n), # suppress default labels
                  col=category_colors[as.character(data_odds$category)],
                  width=0.1,   # slim bars
                  space=0.2,   # reduced space between bars
                  ylim=c(0, max(data_odds$odds_ratio_ci_upper) * 1.3),
                  main="Odds Ratio by Category",
                  xlab="",
                  ylab="Odds Ratio",
                  axes=FALSE)
    axis(2)
    abline(h=0, lwd=1)
    abline(h=1, lty=2, col="red", lwd=2)

    # Add error bars
    arrows(x0=bp,
           y0=data_odds$odds_ratio_ci_lower,
           x1=bp,
           y1=data_odds$odds_ratio_ci_upper,
           angle=90,
           code=3,
           length=0.1)

    # Add x-axis ticks at bar centers
    axis(1, at=bp, labels=FALSE, tck=-0.02)

    # Add 45-degree rotated x-axis labels at bar centers
    text(x=bp, y=par("usr")[3] - 0.05 * max(data_odds$odds_ratio_ci_upper),
         labels=data_odds$category, srt=45, adj=1, xpd=TRUE, cex=1.2)

    # Add value labels with dynamic positioning
    text(x=bp,
         y=data_odds$odds_ratio_ci_upper + 0.1,
         labels=round(data_odds$odds_ratio, 2),
         pos=3, cex=1.5)

    # Add legend
    legend("topleft",
           legend=names(category_colors),
           fill=category_colors,
           cex=1.5,
           bty="n",
           xpd=TRUE)
dev.off()

# 3. Fraction Overlap Plot
# Sort data by percentage
data_perc <- data[order(data$percentage),]
n <- nrow(data_perc)
pdf(file.path(output_dir, "fraction_overlap_plot.pdf"), width=16, height=14)
par(mar=c(10,4,4,6))
bp <- barplot(data_perc$percentage, 
              names.arg=rep('', n), # suppress default labels
              col=category_colors[as.character(data_perc$category)],
              width=0.1,   # slim bars
              space=0.2,   # reduced space between bars
              ylim=c(0, max(data_perc$percentage) * 1.2),
              main="Fraction Overlap by Category",
              xlab="",
              ylab="Percentage (%)",
              axes=FALSE)
axis(2)
abline(h=0, lwd=1)

# Add x-axis ticks at bar centers
axis(1, at=bp, labels=FALSE, tck=-0.02)

# Add 45-degree rotated x-axis labels at bar centers
text(x=bp, y=par("usr")[3] - 0.05 * max(data_perc$percentage),
     labels=data_perc$category, srt=45, adj=1, xpd=TRUE, cex=1.2)

# Calculate a larger dynamic offset (e.g., 5% of the max value)
offset <- 0.05 * max(data_perc$percentage)

# Add value labels with n values, placed just above each bar
text(x=bp,
     y=data_perc$percentage + offset,
     labels=paste0(round(data_perc$percentage, 2), "%\n(n=", data_perc$overlapping_peaks, ")"),
     adj=0.5,
     cex=1.5)

# Add legend
legend("topleft",
       legend=names(category_colors),
       fill=category_colors,
       cex=1.5,
       bty="n",
       xpd=TRUE)
dev.off()




