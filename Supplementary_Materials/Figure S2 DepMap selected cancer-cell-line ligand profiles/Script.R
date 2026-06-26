## =========================================================
## Ligand-level NKG2DL RNA profiles in selected DepMap cell lines
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

## =========================================================
## 1. READ DATA
## =========================================================

df <- read.csv(
  "DepMap_selected_cell_lines_NKG2DL.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## =========================================================
## 2. PREPARE DATA
## =========================================================

gene_order <- c(
  "MICA",
  "MICB",
  "ULBP1",
  "ULBP2",
  "ULBP3",
  "RAET1E",
  "RAET1G",
  "RAET1L"
)

gene_labels <- c(
  "MICA"   = "MICA",
  "MICB"   = "MICB",
  "ULBP1"  = "ULBP1",
  "ULBP2"  = "ULBP2",
  "ULBP3"  = "ULBP3",
  "RAET1E" = "RAET1E / ULBP4",
  "RAET1G" = "RAET1G / ULBP5",
  "RAET1L" = "RAET1L / ULBP6"
)

cell_line_order <- c(
  "HCT 116",
  "A549",
  "PC-3",
  "K-562",
  "SK-N-SH",
  "NCI-H929",
  "Hep G2",
  "MCF7",
  "PANC-1"
)

df_long <- df %>%
  pivot_longer(
    cols = -gene,
    names_to = "cell_line",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    gene_label = gene_labels[as.character(gene)],
    gene_label = factor(gene_label, levels = rev(gene_labels[gene_order])),
    cell_line = factor(cell_line, levels = cell_line_order),
    label = sprintf("%.1f", expression),
    text_color = ifelse(expression >= 3.2, "white", "black")
  )

## =========================================================
## 3. HEATMAP
## =========================================================

p <- ggplot(
  df_long,
  aes(
    x = cell_line,
    y = gene_label,
    fill = expression
  )
) +
  geom_tile(color = NA) +
  geom_text(
    aes(label = label, color = text_color),
    size = 4.0
  ) +
  scale_color_identity() +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#084594",
    limits = c(0, 5.6),
    oob = squish,
    name = "RNA expression value in source file"
  ) +
  labs(
    title = "Ligand-level NKG2DL RNA profiles in selected DepMap cancer cell lines",
    x = "Cell line",
    y = "NKG2D ligand"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 13,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(1.3, "cm"),
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
  filename = "DepMap_selected_cell_lines_NKG2DL_ligand_heatmap.png",
  plot = p,
  width = 12,
  height = 6.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "DepMap_selected_cell_lines_NKG2DL_ligand_heatmap.pdf",
  plot = p,
  width = 12,
  height = 6.2,
  bg = "white"
)

p