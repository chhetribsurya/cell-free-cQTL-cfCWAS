library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

setwd("~/OneDrive - Mass General Brigham/Projects/CWAS/cfChIP")

cf.k27.file <- 'Results/AI_cQTLs/cfChIP.k27.AI.all.sumstats.txt'
wb.k27.file <- 'Results/AI_cQTLs/wb.k27.AI.all.sumstats.txt'

cf.res <- fread(cf.k27.file) %>%
  select(CHR, POS, RSID, NAME, ALL.AF, ALL.BBINOM.P, N.READS) %>%
  filter(RSID != ".")
colnames(cf.res)[5]='ALL.AF.cf'

wb.res <- fread(wb.k27.file) %>%
  select(CHR, POS, RSID, NAME, ALL.AF, ALL.BBINOM.P, N.READS) %>%
  filter(RSID != ".")
colnames(wb.res)[5]='ALL.AF.wb'

df <- inner_join(cf.res, wb.res, by = "RSID") %>%
  select(RSID, CHR = CHR.x, POS = POS.x,
         `cf.Peaks` = NAME.x, ALL.AF.cf, N.READS.x,
         `wb.Peaks` = NAME.y, ALL.AF.wb, N.READS.y)

colnames(df)[6] = 'N.READS.cf'
colnames(df)[9] = 'N.READS.wb'

# Correlation & p-value
ct <- cor.test(df$ALL.AF.cf, df$ALL.AF.wb, method = "pearson", exact = T)
rho_val <- round(unname(ct$estimate), 2)

df <- df %>%
  mutate(group = case_when(
    ALL.AF.cf < 0.5 & ALL.AF.wb < 0.5 ~ "Both allelic fraction < 0.50",
    ALL.AF.cf >= 0.5 & ALL.AF.wb >= 0.5 ~ "Both allelic fraction > 0.50",
    TRUE ~ "Opposite direction"   # one below, one above
  ))

p <- ggplot(df, aes(x = ALL.AF.cf, y = ALL.AF.wb)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 1, color = "gray30") +
  geom_point(aes(color = group), alpha = 0.6) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE) +
  annotate("text", x = 0.05, y = 1, hjust = 0,
           label = paste0("rho == ", rho_val, " * ',' ~~ P < 2.2 %*% 10^-16"),
           parse = TRUE) +
  scale_color_manual(values = c("Both allelic fraction < 0.50" = "red",
                                "Both allelic fraction > 0.50" = "blue",
                                "Opposite direction" = "grey30"),
                     guide = guide_legend(title = NULL)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                     labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                     labels = label_number(accuracy = 0.01)) +
  labs(x = "cfChIP-H3K27ac allele fraction", y = "WBC-H3K27ac allele fraction") +
  theme_classic(base_family = "Helvetica", base_size = 10) +
  theme(
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x  = element_text(size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black")
  )
p
# Save 75 mm x 75 mm vector PDF
ggsave("Manuscripts/Figures/2_e.pdf", plot = p,
       device = "pdf", family = "Helvetica", useDingbats = FALSE,
       width = 130, height = 100, units = "mm")
