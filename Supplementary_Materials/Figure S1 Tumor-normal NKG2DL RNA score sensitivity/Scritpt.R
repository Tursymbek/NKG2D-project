## =========================================================
## NKG2DL tumor-normal contrast sensitivity heatmap
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

## =========================================================
## 1. READ DATA
## =========================================================

df <- read.csv(
  "NKG2DL_tissue_statistics_threshold_sensitivity_0.5_0.6_0.8.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## =========================================================
## 2. PREPARE DATA
## =========================================================

tissue_order <- c(
  "Cervix", "Bladder", "Rectum", "Liver", "Brain",
  "Thyroid", "Kidney", "Colon", "Skin", "Breast",
  "Prostate", "HeadNeck", "Lung", "Uterus", "Ovary"
)

df_plot <- df %>%
  mutate(
    threshold_label = paste0("Purity > ", threshold),
    threshold_label = factor(
      threshold_label,
      levels = rev(c("Purity > 0.5", "Purity > 0.6", "Purity > 0.8"))
    ),
    tissue_simple = factor(tissue_simple, levels = tissue_order),
    stars_clean = ifelse(stars == "ns" | is.na(stars), "", stars),
    label = paste0(sprintf("%.1f", delta_median), stars_clean),
    text_color = ifelse(abs(delta_median) >= 25, "white", "black")
  )

## =========================================================
## 3. HEATMAP
## =========================================================

p <- ggplot(
  df_plot,
  aes(
    x = tissue_simple,
    y = threshold_label,
    fill = delta_median
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = label, color = text_color),
    size = 3.2,
    fontface = "plain"
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-50, 50),
    oob = squish,
    name = "Delta median tumor - normal"
  ) +
  labs(
    title = "Sensitivity of tumor-normal composite NKG2DL RNA contrast to tumor-purity threshold",
    x = "Tissue",
    y = "TCGA tumor-purity filter"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(1.2, "cm"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.line = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

## =========================================================
## 4. SAVE
## =========================================================

ggsave(
  filename = "Figure_S1_NKG2DL_purity_threshold_sensitivity.png",
  plot = p,
  width = 14,
  height = 4.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figure_S1_NKG2DL_purity_threshold_sensitivity.pdf",
  plot = p,
  width = 14,
  height = 4.2,
  bg = "white"
)

p