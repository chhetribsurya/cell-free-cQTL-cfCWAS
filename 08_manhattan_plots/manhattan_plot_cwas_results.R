library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)
library(readr)
# Set working directory (adjust accordingly)
setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

source('../Scripts/functions.R')

## ---- global figure/text settings
base_family <- "Helvetica"
base_size   <- 10           # 10 pt text
w_mm        <- 60           # change to 75 if you want 75x75
h_mm        <- 100           # change to 75 if you want 75x75

theme_set(theme_classic(base_family = base_family, base_size = base_size))

save_pdf <- function(plot, path, w = w_mm, h = h_mm) {
  ggsave(
    filename = path,
    plot = plot,
    device = "pdf",         # vector; no XQuartz needed
    family = base_family,
    useDingbats = FALSE,
    units = "mm",
    width = w, height = h
  )
}

setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

#### Read MPRA file
mpra.df <- read_xlsx('Data/MPRA/df_allelic_activity_cancer_SNV.xlsx')

mpra.df.sub = mpra.df %>%
  select("SNP", "Chr (hg19)", "Position (hg19)", "Ref", "Alt", "Cell types with differential activity", "daSNV") %>%
  filter(!is.na(daSNV))
colnames(mpra.df.sub) = c('RSID', "CHR", "BP", "REF", "ALT", "CELL", "daSNV")

######## Overlap for K4
qtl.df <- read.table('Results/AI_cQTLs/cfChIP.k4.combined.sig.txt',
                     header = TRUE, sep = '\t')
overlap.df <- mpra.df.sub[mpra.df.sub$RSID %in% qtl.df$RSID, ]

cf.names <- read.table('Data/Sample_information/cfChIP_H3K4me3_samples.txt',
                       header = TRUE, sep = "\t",
                       colClasses = c("character", "NULL"))$ID
plot_df <- qtl.df %>% filter(RSID %in% overlap.df$RSID)

SNP <- "rs60143196"

cf.plot.list <- draw_compare_plot(results = plot_df[plot_df$RSID == SNP, ],
                                  names = cf.names,
                                  panel = 'cfChIP-seq',
                                  READ_THRESH = 10,
                                  base_size=10)
p <- cf.plot.list[[1]] + theme(text = element_text(size = base_size, family = base_family, colour = "black"))
q <- cf.plot.list[[2]] + theme(text = element_text(size = base_size, family = base_family, colour = "black"))
p
save_pdf(p, paste0('Manuscripts/Figures/2_c_MPRA.K4.AF.', SNP, '.pdf'))
#save_pdf(q, paste0('Figures/cQTL_summary/MPRA.K4.FOREST.', SNP, '.pdf'))

subtable <- merge(plot_df, mpra.df, by.x="RSID", by.y="SNP") %>%
  select(RSID, CHR, POS, NAME, comb.z, comb.p, Ref, Alt,
         "Cell types with differential activity", "daSNV")
#write_excel_csv(subtable, "Tables/sup5.csv")

######## Overlap for K27
qtl.df <- read.table('Results/AI_cQTLs/cfChIP.k27.combined.sig.txt',
                     header = TRUE, sep = '\t')
overlap.df <- mpra.df.sub[mpra.df.sub$RSID %in% qtl.df$RSID, ]

cf.names <- read.table('Data/Sample_information/cfChIP_H3K27ac_samples.txt',
                       header = TRUE, sep = "\t",
                       colClasses = c("character", "NULL"))$ID
plot_df <- qtl.df %>% filter(RSID %in% overlap.df$RSID)

SNP <- "rs7729529"

cf.plot.list <- draw_compare_plot(results = plot_df[plot_df$RSID == SNP, ],
                                  names = cf.names,
                                  panel = 'cfChIP-seq',
                                  READ_THRESH = 20,
                                  base_size=10)
p <- cf.plot.list[[1]] + theme(text = element_text(size = base_size, family = base_family, colour = "black"))
q <- cf.plot.list[[2]] + theme(text = element_text(size = base_size, family = base_family, colour = "black"))

save_pdf(p, paste0('Manuscripts/Figures/2_d_MPRA.K27.AF.', SNP, '.pdf'))
#save_pdf(q, paste0('Figures/cQTL_summary/MPRA.K27.FOREST.', SNP, '.pdf'))

subtable <- merge(plot_df, mpra.df, by.x="RSID", by.y="SNP") %>%
  select(RSID, CHR, POS, NAME, comb.z, comb.p, Ref, Alt,
         "Cell types with differential activity", "daSNV")
#write_excel_csv(subtable, "Tables/sup6.csv")
