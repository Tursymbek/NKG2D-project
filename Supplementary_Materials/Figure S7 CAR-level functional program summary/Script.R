## =========================================================
## CAR-level functional program summary
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

## =========================================================
## 1. READ DATA
## =========================================================

data_dir <- "."   # если файл в другой папке, укажи путь сюда

df <- read.csv(
  file.path(data_dir, "CAR_program_summary_for_Figure_S7.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## =========================================================
## 2. SETTINGS
## =========================================================

car_order <- c(
  "CD357-CD79B",
  "CD4-FCGR2A",
  "TIM1-FCER1G",
  "HVEM-CD79B",
  "4-1BB-CD79A",
  "Unstim CD28-CD3Z",
  "TIM1-CD3G",
  "CD4-K1",
  "CD223-FCER1G",
  "CD30-CD79B",
  "TIM1-CD79B",
  "4-1BB-FCER1G",
  "CD223-CD79B"
)

program_order <- c(
  "Cytotoxicity",
  "Activation",
  "Proliferation",
  "NKG2DL score"
)

## =========================================================
## 3. PREPARE DATA
## =========================================================

df_long <- df %>%
  pivot_longer(
    cols = all_of(program_order),
    names_to = "Program",
    values_to = "Score"
  ) %>%
  mutate(
    CAR_clean = factor(CAR_clean, levels = rev(car_order)),
    Program = factor(Program, levels = program_order),
    label = sprintf("%.2f", Score),
    text_color = ifelse(abs(Score) >= 1.5, "white", "black")
  )

max_abs <- max(abs(df_long$Score), na.rm = TRUE)

## =========================================================
## 4. PLOT
## =========================================================

p <- ggplot(
  df_long,
  aes(
    x = Program,
    y = CAR_clean,
    fill = Score
  )
) +
  geom_tile(color = NA) +
  geom_text(
    aes(label = label, color = text_color),
    size = 3.6
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = squish,
    name = "Scaled mean program score"
  ) +
  labs(
    title = "CAR-level summary of functional gene programs and NKG2DL score",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 12
    ),
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(1.4, "cm"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.line = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

## =========================================================
## 5. SAVE
## =========================================================

ggsave(
  filename = "Figure_S7_CAR_functional_program_summary.png",
  plot = p,
  width = 9.5,
  height = 7.0,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figure_S7_CAR_functional_program_summary.pdf",
  plot = p,
  width = 9.5,
  height = 7.0,
  bg = "white"
)

p