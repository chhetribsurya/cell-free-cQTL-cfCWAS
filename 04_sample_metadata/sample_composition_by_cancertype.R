setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

library(readxl)
library(dplyr)
library(ggplot2)
library(extrafont)  # To manage fonts if needed
library(tidyr)

theme_set(
  theme_classic(base_family = "Helvetica", base_size = 10)
)

meta.data = read_excel('~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP/Data/Sample_information/meta_data.xlsx') 
k27.cfCWAS.sample = read.table('~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP/Data/Sample_information/cfChIP_H3K27ac_samples.txt', header = T)
k4.cfCWAS.sample = read.table('~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP/Data/Sample_information/cfChIP_H3K4me3_samples.txt', header = T)
cfCWAS.sample = rbind.data.frame(k27.cfCWAS.sample,
                                 k4.cfCWAS.sample
                                 )

meta.data = merge(meta.data, cfCWAS.sample, by.y='ID', by.x = 'study_name')

df_plot = meta.data %>%
  group_by(cancer_type, antibody) %>%
  count()


df_plot = df_plot %>%
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
      TRUE ~ cancer_type  
    )
  )

lev = rev(df_plot[df_plot$antibody=="H3K4me3", ][order(df_plot[df_plot$antibody=="H3K4me3", ]$n),]$cancer_type)
df_plot$cancer_type = factor(df_plot$cancer_type, levels = lev)
df_plot$antibody = factor(df_plot$antibody, levels = c('H3K4me3',
                                                       'H3K27ac'),
                          labels = c('H3K4me3 (N=448)',
                                     'H3K27ac (N=303)'))

colnames(df_plot)[2] ="Assays"

df_plot <- df_plot  %>%
  ungroup() %>%
  complete(cancer_type, Assays, fill = list(n = 0))

# Create the plot
p <- ggplot(df_plot, aes(x = cancer_type, y = n, fill = Assays, colour = Assays, group = Assays, label = n)) +
  geom_col(
    position = position_dodge2(width = 0.9, preserve = "single", padding = 0.2),  # better for zeros
    width = 0.8
  )+
  scale_fill_manual(values = c(
    "H3K4me3 (N=448)" = "#50ad9f",
    "H3K27ac (N=303)" = "#e9c716"
  ))+
  scale_color_manual(values = c(
    "H3K4me3 (N=448)" = "#50ad9f",
    "H3K27ac (N=303)" = "#e9c716"
  ))+
  scale_y_continuous(expand = c(0, 0))+
  theme(
    axis.text.x  = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    legend.text  = element_text(size = 10, colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  xlab("Cancer Types") +
  ylab("Sample size")

p

# Save the plot as a vector PDF file
ggsave(filename = "Manuscripts/Figures/1_a.pdf",
       plot = p,
       device = pdf,  # Ensure vector output
       units = "mm",
       family = "Helvetica",      # use mac Helvetica
       useDingbats = FALSE,       # avoids weird font substitutions
       width = 120, height = 100)




