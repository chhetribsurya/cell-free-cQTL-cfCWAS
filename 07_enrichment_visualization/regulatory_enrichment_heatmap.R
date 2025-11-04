library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

#--------- Helpers ---------#

split_num <- function(x) {
  if (is.na(x) || x == "") return(numeric(0))
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

# Build a long dataframe for one RSID (safe if one of the conditions is missing)
build_plotmat <- function(df, df_names, rsid) {
  i <- which(df$RSID == rsid)
  if (length(i) != 1) stop("RSID not found or not unique: ", rsid)
  
  n0 <- length(split_num(df$IND.C0[i]))
  n1 <- length(split_num(df$IND.C1[i]))
  
  df0 <- if (n0 > 0) {
    tibble(
      rsid      = df$RSID[i],
      condition = "WBC ChIP",
      ind       = split_num(df$IND.C0[i]),
      ref       = split_num(df$IND.C0.COUNT.REF[i]),
      alt       = split_num(df$IND.C0.COUNT.ALT[i])
    )
  } else NULL
  
  df1 <- if (n1 > 0) {
    tibble(
      rsid      = df$RSID[i],
      condition = "cfChIP",
      ind       = split_num(df$IND.C1[i]),
      ref       = split_num(df$IND.C1.COUNT.REF[i]),
      alt       = split_num(df$IND.C1.COUNT.ALT[i])
    )
  } else NULL
  
  plotmat <- bind_rows(df0, df1)
  
  # Map individual indices to sample IDs; keep a stable order
  plotmat <- plotmat %>%
    mutate(ind = factor(df_names[ind], levels = df_names[order(df_names)])) %>%
    pivot_longer(c(ref, alt), names_to = "allele", values_to = "counts") %>%
    mutate(
      allele    = factor(str_to_upper(allele), levels = c("REF","ALT")),
      condition = factor(condition, levels = c("WBC ChIP","cfChIP"))
    )
  
  plotmat
}

plot_ai_for_rsid <- function(df, df_names, rsid) {
  plotmat <- build_plotmat(df, df_names, rsid)
  i <- which(df$RSID == rsid)
  if (length(i) != 1) stop("RSID not found or not unique: ", rsid)
  
  # Values for annotations
  cf.AF <- round(df$C1.AF[i], 2)
  wb.AF <- round(df$C0.AF[i], 2)
  p.cf  <- df$C1.BBINOM.P[i]
  p.wb  <- df$C0.BBINOM.P[i]
  
  # Format p-values for plotmath
  p_to_plotmath_rhs <- function(p, thresh = .Machine$double.eps) {
    if (is.na(p)) return("NA")
    if (p > 0.05) {
      # Show as rounded number only
      return(sprintf("%.2f", p))
    } else if (p < thresh) {
      # Very small p-values
      return("2.2 %*% 10^-16")  # like cor.test output
    } else {
      # Scientific notation
      e <- floor(log10(p))
      m <- p / (10^e)
      paste0(sprintf("%.2f", m), " %*% 10^", e)
    }
  }
  
  # Labels using atop (two lines: AF and p)
  cf_label <- paste0(
    "atop(AF[cf] == ", sprintf("%.2f", cf.AF), ", p[cf] == ", p_to_plotmath_rhs(p.cf), ")"
  )
  wb_label <- paste0(
    "atop(AF[wb] == ", sprintf("%.2f", wb.AF), ", p[wb] == ", p_to_plotmath_rhs(p.wb), ")"
  )
  
  # Totals per bar for placement
  totals <- plotmat %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(total = sum(counts), .groups = "drop")
  
  # Offset for labels
  offset <- 0.04 * max(totals$total, na.rm = TRUE)
  
  ann <- tibble::tibble(
    condition = factor(c("cfChIP", "WBC ChIP"), levels = c("WBC ChIP", "cfChIP")),
    y = c(
      totals$total[match("cfChIP", totals$condition)],
      totals$total[match("WBC ChIP", totals$condition)]
    ) + offset,
    label = c(cf_label, wb_label)
  )
  capitalize_first <- function(x) {
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }
  ggplot(plotmat, aes(x = condition, y = counts, fill = allele, colour = allele)) +
    geom_col(position = "stack", width = 0.8) +
    geom_text(
      data = ann,
      aes(x = condition, y = y, label = label),
      inherit.aes = FALSE,
      parse = TRUE,       # necessary for atop()
      vjust = 0,
      size = 8 / .pt,
      colour = "black"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = c("REF" = "#1f78b4", "ALT" = "#e31a1c"), name = "Allele") +
    scale_color_manual(values = c("REF" = "#1f78b4", "ALT" = "#e31a1c"), guide = "none") +
    labs(
      title = capitalize_first(df$NAME[i]),
      x = NULL, y = "Read counts"
    ) +
    theme_classic(base_family = "Helvetica", base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, colour = "black"),
      legend.text  = element_text(size = 10, colour = "black"),
      axis.text    = element_text(size = 10, colour = "black"),
      axis.title   = element_text(size = 10, colour = "black"),
      plot.title   = element_text(size = 10, colour = "black")
    )
}




# Compute summed counts per condition for multiple RSIDs (tidy, no nested if/else)
sum_counts_by_condition <- function(df, df_names, rsids) {
  rsids %>%
    map_df(function(snp) {
      pm <- build_plotmat(df, df_names, snp)
      pm %>%
        group_by(rsid, condition) %>%
        summarise(total = sum(counts), .groups = "drop") %>%
        pivot_wider(names_from = condition, values_from = total,
                    names_prefix = "total.") %>%
        mutate(rsid = snp)
    })
}

#--------- Load and prepare data ---------#

df <- fread("Results/cf_wb_AI/cf_peaks/results.all.txt")
df_names <- fread("Results/cf_wb_AI/samples.txt")$ID  # 1-based indices in your data

#--------- Example: single RSID plot ---------#
#c:"rs2239677", d:"rs9526983"
rsid <- "rs9526983"
p <- plot_ai_for_rsid(df, df_names, rsid)
print(p)

ggsave(
  filename = file.path("Manuscripts/Figures",
                       paste0("3_c_", rsid, ".context_dependent_peak.", df$NAME[df$RSID == rsid], ".AI.pdf")),
  plot   = p,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  units  = "mm", width = 60, height = 100   # adjust size as needed
)

#--------- Example: batch counts summary ---------#
# Suppose you have a vector of RSIDs in `tmp$RSID`:
# tmp <- data.frame(RSID = c("rs9526983", "rsXXXX", ...))
# counts_summary <- sum_counts_by_condition(df, df_names, tmp$RSID)
# head(counts_summary)
