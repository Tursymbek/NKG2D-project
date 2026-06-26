## =========================================================
## ENCODE promoter-overlap summary for candidate NKG2DL regulators
## =========================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})

## =========================================================
## 1. PATHS
## =========================================================

data_dir <- "."   # поменяй на папку с файлами, если нужно

matrix_file <- file.path(
    data_dir,
    "MATRIX_TF_x_NKG2DL_total_peaks.csv"
)

## =========================================================
## 2. READ DATA
## =========================================================

mat <- read.csv(
    matrix_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

## =========================================================
## 3. PREPARE DATA
## =========================================================

ligand_order <- c(
    "MICA",
    "MICB",
    "ULBP1",
    "ULBP2",
    "ULBP3",
    "RAET1E",
    "RAET1G",
    "RAET1L"
)

ligand_labels <- c(
    "MICA"   = "MICA",
    "MICB"   = "MICB",
    "ULBP1"  = "ULBP1",
    "ULBP2"  = "ULBP2",
    "ULBP3"  = "ULBP3",
    "RAET1E" = "RAET1E / ULBP4",
    "RAET1G" = "RAET1G / ULBP5",
    "RAET1L" = "RAET1L / ULBP6"
)

mat <- mat %>%
    select(TF, all_of(ligand_order))

tf_summary <- mat %>%
    mutate(total_overlapping_promoter_peaks = rowSums(across(all_of(ligand_order)), na.rm = TRUE)) %>%
    arrange(desc(total_overlapping_promoter_peaks)) %>%
    slice_head(n = 20) %>%
    mutate(
        TF = factor(TF, levels = rev(TF))
    )

ligand_summary <- mat %>%
    pivot_longer(
        cols = all_of(ligand_order),
        names_to = "ligand",
        values_to = "n_overlapping_peaks"
    ) %>%
    group_by(ligand) %>%
    summarise(
        n_candidate_TFs_with_overlap = sum(n_overlapping_peaks > 0, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        ligand_label = ligand_labels[ligand],
        ligand_label = factor(
            ligand_label,
            levels = rev(ligand_labels[ligand_order])
        )
    )

## =========================================================
## 4. PANEL A
## =========================================================

p1 <- ggplot(
    tf_summary,
    aes(
        x = total_overlapping_promoter_peaks,
        y = TF
    )
) +
    geom_col(
        fill = "#9B78C8",
        width = 0.75
    ) +
    scale_x_continuous(
        limits = c(0, 960),
        breaks = seq(0, 900, 200),
        expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
        title = "A. Candidate TFs with largest total promoter overlap",
        x = "Total overlapping promoter peaks",
        y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 16,
            hjust = 0
        ),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

## =========================================================
## 5. PANEL B
## =========================================================

p2 <- ggplot(
    ligand_summary,
    aes(
        x = n_candidate_TFs_with_overlap,
        y = ligand_label
    )
) +
    geom_col(
        fill = "#F39A32",
        width = 0.75
    ) +
    scale_x_continuous(
        limits = c(0, 52),
        breaks = seq(0, 50, 10),
        expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
        title = "B. TF coverage by NKG2DL promoter",
        x = "Candidate TFs with at least one overlap",
        y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 16,
            hjust = 0
        ),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 13),
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

## =========================================================
## 6. COMBINE
## =========================================================

p_final <- (p1 | p2) +
    plot_layout(widths = c(1.25, 1)) +
    plot_annotation(
        title = "ENCODE promoter-overlap summary for candidate NKG2DL regulators",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 20,
                hjust = 0.5
            )
        )
    )

## =========================================================
## 7. SAVE
## =========================================================

ggsave(
    filename = "Figure_S4_ENCODE_promoter_overlap_summary.png",
    plot = p_final,
    width = 13.5,
    height = 6.3,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "Figure_S4_ENCODE_promoter_overlap_summary.pdf",
    plot = p_final,
    width = 13.5,
    height = 6.3,
    bg = "white"
)

p_final