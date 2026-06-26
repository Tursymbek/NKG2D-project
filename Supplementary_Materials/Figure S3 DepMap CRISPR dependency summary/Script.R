## =========================================================
## DepMap CRISPR dependency summary for NKG2D ligand genes
## =========================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})

## =========================================================
## 1. READ DATA
## =========================================================

df <- read.delim(
    "DepMap_NKG2DL_CRISPRGeneEffect_summary.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

## =========================================================
## 2. PREPARE DATA
## =========================================================

gene_order <- c(
    "MICA", "MICB", "ULBP1", "ULBP2",
    "ULBP3", "ULBP4", "ULBP5", "ULBP6"
)

df_plot <- df %>%
    mutate(
        Gene = factor(Gene, levels = rev(gene_order))
    )

## =========================================================
## 3. PANEL A: MEDIAN CHRONOS GENE-EFFECT SCORE
## =========================================================

p1 <- ggplot(
    df_plot,
    aes(
        x = median_effect,
        y = Gene
    )
) +
    geom_col(
        fill = "#6F8FB9",
        width = 0.8
    ) +
    geom_vline(
        xintercept = 0,
        color = "black",
        linewidth = 0.5
    ) +
    geom_vline(
        xintercept = -0.5,
        color = "#E64B35",
        linetype = "dashed",
        linewidth = 0.7
    ) +
    scale_x_continuous(
        limits = c(-0.53, 0.19),
        breaks = seq(-0.5, 0.1, 0.1),
        expand = expansion(mult = c(0.01, 0.03))
    ) +
    labs(
        title = "A. Median effect near neutral",
        x = "Median Chronos gene-effect score",
        y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 16,
            hjust = 0
        ),
        axis.text.x = element_text(
            size = 12,
            color = "black"
        ),
        axis.text.y = element_text(
            size = 12,
            color = "black"
        ),
        axis.title.x = element_text(size = 13),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.grid.major.x = element_line(
            color = "grey88",
            linewidth = 0.4
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

## =========================================================
## 4. PANEL B: PERCENT BELOW -0.5
## =========================================================

p2 <- ggplot(
    df_plot,
    aes(
        x = percent_below_minus_0_5,
        y = Gene
    )
) +
    geom_col(
        fill = "#7EAD5F",
        width = 0.8
    ) +
    scale_x_continuous(
        limits = c(0, 2.52),
        breaks = seq(0, 2.5, 0.5),
        expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
        title = "B. Rare partial-dependency calls",
        x = "Cell lines below -0.5 (%)",
        y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 16,
            hjust = 0
        ),
        axis.text.x = element_text(
            size = 12,
            color = "black"
        ),
        axis.text.y = element_text(
            size = 12,
            color = "black"
        ),
        axis.title.x = element_text(size = 13),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.grid.major.x = element_line(
            color = "grey88",
            linewidth = 0.4
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

## =========================================================
## 5. COMBINE PANELS
## =========================================================

p_final <- (p1 | p2) +
    plot_annotation(
        title = "DepMap CRISPR dependency summary for NKG2D ligand genes",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 20,
                hjust = 0.5
            )
        )
    )

## =========================================================
## 6. SAVE
## =========================================================

ggsave(
    filename = "DepMap_NKG2DL_CRISPR_dependency_summary.png",
    plot = p_final,
    width = 13,
    height = 6,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "DepMap_NKG2DL_CRISPR_dependency_summary.pdf",
    plot = p_final,
    width = 13,
    height = 6,
    bg = "white"
)

p_final