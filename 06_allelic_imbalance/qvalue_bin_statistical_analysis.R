setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")
library(qvalue)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Hmisc)
library(VennDiagram)
library(readr)
library(data.table)
library(rstatix)
library(ggpubr)
library(forcats)
source("../Scripts/functions.R")
set.seed(123)
theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

SIG=0.05 # Q value significance threshold
READ_THRESH=20 # SNPs must have at least this many reads used in AI calculation to be tested for significance
BAL_THRESH=100 # SNPs must have at least this many reads used in AI calculation to be added to list of covered but balanced SNPs


##############
all.results=fread("Results/cf_wb_AI/cf_peaks/results.all.txt",sep="\t",header=TRUE,colClasses=c("character", "numeric", "character", "numeric", "numeric", "character", rep("numeric",10), rep("character",6)))

IND.C0=get_counts(all.results$IND.C0)
IND.C0.COUNT.REF=get_counts(all.results$IND.C0.COUNT.REF)
all.results$IND.C0.COUNT.REF = unlist(lapply(IND.C0.COUNT.REF, FUN = sum))

IND.C0.COUNT.ALT=get_counts(all.results$IND.C0.COUNT.ALT)
all.results$IND.C0.COUNT.ALT = unlist(lapply(IND.C0.COUNT.ALT, FUN = sum))

all.results$C0.COUNT = all.results$IND.C0.COUNT.REF + all.results$IND.C0.COUNT.ALT

IND.C1=get_counts(all.results$IND.C1)
IND.C1.COUNT.REF=get_counts(all.results$IND.C1.COUNT.REF)
all.results$IND.C1.COUNT.REF = unlist(lapply(IND.C1.COUNT.REF, FUN = sum))

IND.C1.COUNT.ALT=get_counts(all.results$IND.C1.COUNT.ALT)
all.results$IND.C1.COUNT.ALT = unlist(lapply(IND.C1.COUNT.ALT, FUN = sum))

all.results$C1.COUNT = all.results$IND.C1.COUNT.REF + all.results$IND.C1.COUNT.ALT

all.results$count.ratio = all.results$C1.COUNT/(all.results$C0.COUNT)


results=all.results[(all.results$C1.COUNT >= READ_THRESH|all.results$C0.COUNT >= READ_THRESH),]
results$POSm1=results$POS-1

# get qvalues values
results$BBINOM.Q=qvalue(results$ALL.BBINOM.P)$qvalue
results$BBINOM.C0.Q=qvalue(results$C0.BBINOM.P)$qvalue
results$BBINOM.C1.Q=qvalue(results$C1.BBINOM.P)$qvalue
results$DIFF.Q=qvalue(results$DIFF.BBINOM.P)$qvalue

#####################
#####################
sig.C1.specific <- subset(
  results,
  (BBINOM.C1.Q < SIG & BBINOM.C0.Q >= SIG) | (BBINOM.C1.Q >= SIG & BBINOM.C0.Q >= SIG)
)

# Helpers that return plotmath-ready CHARACTER strings
sci <- function(x) paste0("1 %*% 10^", x)          # e.g. "1 %*% 10^-5"
lt  <- function(x) paste0("p < ", x)               # "p < 1 × 10^-10"
rng <- function(a, b) paste0(a, " <= p ~ '<' ~ ", b)  # "a <= p < b" (middle '<' as text)

# Assign plotmath-style bins
df <- sig.C1.specific %>%
  mutate(
    p_bin = case_when(
      BBINOM.C1.Q < 1e-10                         ~ lt(sci("-10")),
      BBINOM.C1.Q >= 1e-10 & BBINOM.C1.Q < 1e-5   ~ rng(sci("-10"), sci("-5")),
      BBINOM.C1.Q >= 1e-5  & BBINOM.C1.Q < 1e-2   ~ rng(sci("-5"),  sci("-2")),
      BBINOM.C1.Q >= 1e-2  & BBINOM.C1.Q < 0.05   ~ rng(sci("-2"),  "0.05"),
      BBINOM.C1.Q >= 0.05                         ~ "p >= 0.05",
      TRUE ~ NA_character_
    ),
    p_bin = factor(
      p_bin,
      levels = c(
        lt(sci("-10")),
        rng(sci("-10"), sci("-5")),
        rng(sci("-5"),  sci("-2")),
        rng(sci("-2"),  "0.05"),
        "p >= 0.05"
      )
    )
  )


# ANOVA on bins except the non-significant bin
anova_result <- aov(log10(count.ratio) ~ p_bin, data = subset(df, p_bin != "p >= 0.05"))
summary(anova_result)

#Non-parametric test
kw <- kruskal_test(subset(df, p_bin != "p >= 0.05"), log10(count.ratio) ~ p_bin)
kw

# Subsample 20% only for the non-significant bin
set.seed(123)
df <- df %>%
  group_by(p_bin) %>%
  mutate(sampled = if_else(p_bin == "p >= 0.05", row_number() %% 50 == 0, TRUE)) %>%
  ungroup() %>%
  filter(sampled)

# ---- Boxplot ----
df_pw <- df %>%
  filter(#p_bin != "p >= 0.05",
    is.finite(count.ratio),
    count.ratio > 0) %>%
  mutate(
    log10ratio = log10(count.ratio),
    p_bin = fct_drop(p_bin)
  )

lv <- levels(df_pw$p_bin)
comparisons <- list(
  c("p < 1 %*% 10^-10","1 %*% 10^-10 <= p ~ '<' ~ 1 %*% 10^-5"),
  c("1 %*% 10^-10 <= p ~ '<' ~ 1 %*% 10^-5","1 %*% 10^-5 <= p ~ '<' ~ 1 %*% 10^-2"),
  c("1 %*% 10^-5 <= p ~ '<' ~ 1 %*% 10^-2", "1 %*% 10^-2 <= p ~ '<' ~ 0.05"),
  c("1 %*% 10^-2 <= p ~ '<' ~ 0.05", "p >= 0.05")
)

p.count.boxplot <- ggplot(df, aes(x = p_bin, y = log10(count.ratio), fill = p_bin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_jitter(position = position_jitter(width = 0.35, height = 0.10), alpha = 0.25, size = 0.5) +
  labs(
    x = "P-value bin (Non-WBC cQTLs)",
    y = "log10(cfChIP / WBC read-count ratio)"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  scale_fill_brewer(palette = "Blues", direction = -1, name = NULL) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  coord_cartesian(ylim = c(-2.7, 2.7), clip = "off") +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title   = element_text(size = 10, colour = "black"),
    plot.margin  = margin(5, 5, 18, 5)
  )+
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     p.adjust.method = "BH",
                     label = "p.signif",
                     step.increase = 0.03,
                     tip.length = 0.01)

print(p.count.boxplot)

ggsave('Manuscripts/Figures/ext_2_b.pdf', p.count.boxplot,
       width = 140, height = 120, units = 'mm',
       device = "pdf", family = "Helvetica", useDingbats = FALSE)
