library(qvalue)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Hmisc)
library(VennDiagram)
library(readr)
library(data.table)
setwd("~/Projects/cfChIP")
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
results$log_ratio = log10(results$count.ratio)
#####################

# cfChIP-significant & WBC-non-significant
sig.C1.specific <- subset(results, BBINOM.C1.Q < SIG & BBINOM.C0.Q >= SIG)

# Prepare x data safely
sig.C1.specific <- sig.C1.specific |>
  mutate(log_ratio = log10(count.ratio)) |>
  filter(is.finite(log_ratio))

# Choose sensible x-limits (2.5th–97.5th pct) to avoid long tails
xr <- quantile(sig.C1.specific$log_ratio, c(0.025, 0.975), na.rm = TRUE)

p.count.dist <- ggplot(sig.C1.specific, aes(x = log_ratio)) +
  geom_histogram(binwidth = 0.1, fill = "#1f78b4", color = "#1f78b4") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
  geom_vline(xintercept = mean(sig.C1.specific$log_ratio), linetype = 1, linewidth = 0.6, colour = "firebrick") +
  labs(
    x = "log10(cfChIP / WBC read-count ratio)",
    y = "Count",
    title = NULL
  ) +
  coord_cartesian(xlim = c(-2,2)) +
  theme(
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    #plot.margin = margin(1,1,1,1)
  )

p.count.dist

# Save as vector PDF
ggsave(
  filename = "Manuscripts/Figures/ext_2a.pdf",
  plot = p.count.dist,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  width = 90, height = 100, units = "mm"
)


###################
###Histogram for all the overlapped peaks
# Prepare x data safely
results <- results |>
  mutate(log_ratio = log10(count.ratio)) |>
  filter(is.finite(log_ratio))

# Choose sensible x-limits (2.5th–97.5th pct) to avoid long tails
xr <- quantile(results$log_ratio, c(0.025, 0.975), na.rm = TRUE)

p.count.dist.all <- ggplot(results, aes(x = log_ratio)) +
  geom_histogram(binwidth = 0.1, fill = "#1f78b4", color = "#1f78b4") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
  geom_vline(xintercept = mean(results$log_ratio), linetype = 1, linewidth = 0.6, colour = "firebrick") +
  labs(
    x = "log10(cfChIP / WBC read-count ratio)",
    y = "Count",
    title = "All tested peaks across cfChIP and WBC"
  ) +
  coord_cartesian(xlim = c(-2,2)) +
  theme(
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    #plot.margin = margin(1,1,1,1)
  )

p.count.dist.all
# Save as vector PDF
ggsave(
  filename = "Figures/cf_wb_AI/p.count.dist.all.pdf",
  plot = p.count.dist.all,
  device = "pdf", family = "Helvetica", useDingbats = FALSE,
  width = 100, height = 100, units = "mm"
)


##################summary
diff.df = results %>%
  #filter(RSID %in% tmp$RSID| RSID %in% tmp2$RSID ) %>%
  filter( BBINOM.C0.Q>0.05, BBINOM.C1.Q<0.05)
#%>%
#  filter(!is.na(DIFF.Q), log_ratio>=log10(1/1.5), log_ratio<=log10(1.5)) 

dim(diff.df)


diff.df = results %>%
  #filter(RSID %in% tmp$RSID | RSID %in% tmp2$RSID ) %>%
  filter( BBINOM.C0.Q>0.05, BBINOM.C1.Q<0.05) %>%
  filter(DIFF.Q< 0.05)

#filter(!is.na(DIFF.Q), log_ratio>=log10(1/1.5), log_ratio<=log10(1.5)) %>%
#  filter( BBINOM.C0.Q>0.05, BBINOM.C1.Q<0.05)

dim(diff.df)

diff.df = results %>%
  filter(!is.na(DIFF.Q), log_ratio>=log10(1/1.5), log_ratio<=log10(1.5)) %>%
  filter(DIFF.Q<0.05)

dim(diff.df)

sum(diff.df$DIFF.Q<0.05)

sum(diff.df$DIFF.Q<0.05)/dim(diff.df)[1]

sum(diff.df$DIFF.BBINOM.P<0.05)

sum(diff.df$DIFF.BBINOM.P<0.05)/dim(diff.df)[1]


tmp = fread('Results/AI_cQTLs/cfChIP.k27.combined.sig.txt') %>%
  mutate(q=qvalue(pval)$qvalue) %>%
  filter(q<0.05)

tmp2 = fread('Results/AI_cQTLs/WBC.k27.combined.sig.txt') %>%
  mutate(q=qvalue(pval)$qvalue) %>%
  filter(q>0.05)
