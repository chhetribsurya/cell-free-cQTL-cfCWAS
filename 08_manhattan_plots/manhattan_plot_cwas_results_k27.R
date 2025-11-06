library(data.table)
library(dplyr)
library(ggplot2)

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

setwd("~/Projects/cfChIP")

cf.k27.file <- 'Results/AI_cQTLs/cfChIP.k27.cQTLs.sumstats.txt'
wb.k27.file <- 'Results/AI_cQTLs/wb.k27.cQTLs.sumstats.txt'

cf.res <- fread(cf.k27.file) %>% mutate(SNP = paste0(chr, ":", snp.pos))
wb.res <- fread(wb.k27.file) %>% mutate(SNP = paste0(chr, ":", snp.pos))

df <- inner_join(cf.res, wb.res, by = "SNP") %>%
  transmute(
    SNP.POS = SNP,
    CHR = chr.x,
    `cf.Peaks.start` = start.x, `cf.Peaks.end` = end.x, beta.cf = beta.x,
    `wb.Peaks.start` = start.y, `wb.Peaks.end` = end.y, beta.wb = beta.y
  )

df <- df %>%
  mutate(group = case_when(
    beta.cf < 0 & beta.wb < 0 ~ "Both effect sizes < 0",
    beta.cf >= 0 & beta.wb >= 0 ~ "Both effect sizes > 0",
    TRUE ~ "Opposite direction"   # one below, one above
  ))

ct <- cor.test(df$beta.cf, df$beta.wb, method = "pearson")
rho_val <- round(unname(ct$estimate), 2)

# Plot with superscripted P using plotmath
p <- ggplot(df, aes(beta.cf, beta.wb)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 1, color = "gray30") +
  geom_point(aes(color = group), alpha = 0.6) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE, linewidth = 1) +
  labs(x = "cfChIP-H3K27ac cQTL effect size", y = "WBC-H3K27ac cQTL effect size") +
  annotate("text", x = -1, y = max(df$beta.wb, na.rm = TRUE),
           hjust = 0, vjust = 1,
           # If extremely small, print "< 2.2 × 10^-16"; else show scientific with superscript
           label = if (ct$p.value < .Machine$double.eps) {
             paste0("rho == ", rho_val, " * ',' ~~ P < 2.2 %*% 10^-16")
           } else {
             # convert ct$p.value to mantissa × 10^exp for superscript
             {
               e <- floor(log10(ct$p.value))
               m <- ct$p.value / (10^e)
               paste0("rho == ", rho_val, " * ',' ~~ P == ",
                      sprintf("%.2f", m), " %*% 10^", e)
             }
           },
           parse = TRUE,
           #size = 10/.pt,
           colour = "black") +
  scale_color_manual(values = c("Both effect sizes < 0" = "red",
                                "Both effect sizes > 0" = "blue",
                                "Opposite direction" = "grey30"),
                     guide = guide_legend(title = NULL)) +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-2, 2, by = 0.5),
                     labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-1.6, 1.9), breaks = seq(-2, 2, by = 0.5),
                     labels = label_number(accuracy = 0.01)) +
  theme_classic(base_family = "Helvetica", base_size = 10)+
  theme(
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x  = element_text(size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black")
  )


print(p)

ggsave("Manuscripts/Figures/2_f.pdf", p,
       device = "pdf", family = "Helvetica", useDingbats = FALSE,
       width = 130, height = 100, units = "mm")
