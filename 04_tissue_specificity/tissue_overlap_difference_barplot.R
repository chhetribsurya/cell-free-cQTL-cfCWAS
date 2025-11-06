library(dplyr)
library(ggplot2)
library(forcats)
library(readr)     # for write_excel_csv
library(extrafont) # if you want Helvetica in PDF

theme_set(theme_classic(base_family = "Helvetica", base_size = 10))

setwd("~/Projects/cfChIP")

state <- 'EnhA1_EnhA2_TssA_EnhG1_EnhG2'
result.file <- file.path(
  '~/Projects/cfChIP/Results/EpiMap_overlap',
  paste0(state, '_combined_results.csv')
)

df <- read.csv(result.file, header = TRUE) %>%
  mutate(
    CF_freq  = CF_count / Total_count,
    WB_freq  = WB_count / Total_count,
    ratio    = CF_count / WB_count,
    delta    = (CF_freq - WB_freq) * 100,
    logratio = log(CF_count / WB_count)
  )

# Output table sup1
out.cf <- df %>%
  select(EID, CF_count, Total_count, NAME, CF_freq) %>%
  rename(
    EID_overlap_count   = CF_count,
    EID_total_count     = Total_count,
    EID_NAME            = NAME,
    PERCENTAGE          = CF_freq
  )
write_excel_csv(out.cf, "Tables/sup1.csv")

# Output table sup3
out.ratio <- df %>%
  arrange(ratio) %>%
  rename(
    EID_overlap_count_WBC          = WB_count,
    EID_overlap_count_cfChIP       = CF_count,
    EID_total_count                = Total_count,
    EID_NAME                       = NAME,
    EID_overlap_percentage_cfChIP  = CF_freq,
    EID_overlap_percentage_WBC     = WB_freq,
    Ratio_cfChIP_to_WBC            = ratio,
    DIFF_cfChIP_to_WBC             = delta,
    Logratio_cfChIP_to_WBC         = logratio
  )
write_excel_csv(out.ratio, "Tables/sup3.csv")

# Prepare df.ratio for plot
df.ratio <- bind_rows(
  df %>% arrange(ratio) %>% slice(1:15),
  df %>% arrange(desc(ratio)) %>% slice(1:10)
) %>%
  mutate(
    sign  = if_else(delta > 0, "cfChIP > WBC ChIP", "WBC ChIP > cfChIP"),
    NAME  = fct_reorder(NAME, delta)
  )

# Plot
p <- ggplot(df.ratio, aes(x = NAME, y = delta, fill = sign)) +
  geom_col(width = 0.8, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") +
  scale_fill_manual(
    values = c("cfChIP > WBC ChIP" = "#1f78b4",
               "WBC ChIP > cfChIP" = "#e31a1c"),
    guide = guide_legend(position = "bottom"),
    name  = NULL
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(y = "Difference in % overlap", x = "Tissues") +
  theme(
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 5), colour = "black"),
    axis.text.x  = element_text(size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    legend.text = element_text(size = 9, colour = "black"),
    legend.title = element_blank(),
    plot.margin = margin(t = 0, r = 100, b = 0, l = 0, unit = "pt")
    )
p
# Save as vector PDF with embedded fonts
ggsave(
  filename = paste0('Manuscripts/Figures/1_c_', state, '.pdf'),
  plot = p,
  device = "pdf",            # base quartz PDF on macOS
  family = "Helvetica",      # use mac Helvetica
  useDingbats = FALSE,       # avoids weird font substitutions
  units = "mm",
  width = 150, height = 120
)

