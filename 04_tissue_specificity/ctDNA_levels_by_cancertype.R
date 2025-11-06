setwd("~/Projects/cfChIP")

library(dplyr)
library(ggplot2)
library(readxl)
library(extrafont)  # For managing fonts
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 10)
)
# Import your data
meta.data <- read_excel('Data/Sample_information/meta_data.xlsx')
k27.cfCWAS.sample <- read.table('Data/Sample_information/cfChIP_H3K27ac_samples.txt', header = TRUE)
k4.cfCWAS.sample <- read.table('Data/Sample_information/cfChIP_H3K4me3_samples.txt', header = TRUE)
cfCWAS.sample <- rbind.data.frame(k27.cfCWAS.sample, k4.cfCWAS.sample)

meta.data <- merge(meta.data, cfCWAS.sample, by.y = 'ID', by.x = 'study_name') %>%
  mutate(ctDNA = as.numeric(ctDNA))

meta.data <- meta.data %>%
  mutate(
    cancer_type = case_when(
      cancer_type == "Breast" ~ "BC",
      cancer_type == "Colorectal" ~ "CRC",
      cancer_type == "Esophageal" ~ "EC",
      cancer_type == "Hepatocellular" ~ "HCC",
      cancer_type == "Merkel cell" ~ "MCC",
      cancer_type == "Non-Small Cell Cancer" ~ "NSCLC",
      cancer_type == "Ovarian" ~ "OC",
      cancer_type == "Prostate" ~ "PRAD",
      cancer_type == "Small Cell Lung Cancer" ~ "SCLC",
      cancer_type == "Small cell bladder" ~ "SCBC",
      cancer_type == "Thymic" ~ "TC",
      cancer_type == "Renal" ~ "RCC",
      cancer_type == "Neuroendocrine Prostate Cancer" ~ "NEPC",
      cancer_type == "Melanoma" ~ "MM",
      cancer_type == "Glioma" ~ "GBM",
      TRUE ~ cancer_type  # Keeps other values unchanged
    )
  )

df.median <- meta.data %>%
  group_by(cancer_type) %>%
  summarize(Median = mean(ctDNA, na.rm = TRUE),
            MAX = max(ctDNA, na.rm = TRUE),
            MIN = min(ctDNA, na.rm = TRUE))

levs <- rev(df.median[order(df.median$Median), ]$cancer_type)
meta.data$cancer_type <- factor(meta.data$cancer_type, levels = levs)

# Create the boxplot with jitter
p <- ggplot(meta.data, aes(y = ctDNA, x = cancer_type)) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "orange") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 0.2) +
  theme_classic(base_family = "Helvetica", base_size = 7) + 
  theme(
    axis.text.x  = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title   = element_text(size = 10, colour = "black"),
    plot.title   = element_text(size = 10, colour = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 7, colour = "black"),
    legend.position = "right"
  ) + 
  labs(fill = "Epigenomic assays") +
  ylab("ctDNA Levels") + 
  xlab("Cancer Types")


p
# Save the figure as a vector PDF file
ggsave("Manuscripts/Figures/1_b.pdf",
       plot = p,
       device = pdf,  # Vector-based output for line art
       width = 120, height = 100, dpi = 300, units = "mm")
