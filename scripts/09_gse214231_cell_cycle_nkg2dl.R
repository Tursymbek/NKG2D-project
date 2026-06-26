# -----------------------------------------------------------------------------
# GSE214231 CAR cell-cycle and NKG2DL pie-scatter analysis
#
# Original archived script: Castellanos-Rueda et al. [76]/C, D - Cell Cycle NKG2DL.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

## =========================================================
## NKG2DL and Cell Cycle
## Includes Unstim CD28-CD3Z as control
## =========================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(scatterpie)
    library(grid)
})


## =========================================================
## 1. LOAD DATA
## =========================================================

seurat_T <- readRDS("seurat_T_with_starCAT.rds")

meta <- read.csv(
    "metadata_cells_with_starCAT.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

scores <- read.table(
    "GSE214231_Tcells.scores.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

## barcode fix
if (colnames(meta)[1] == "") {
    colnames(meta)[1] <- "cell_barcode"
}

if (!"cell_barcode" %in% colnames(meta)) {
    meta$cell_barcode <- as.character(meta[[1]])
} else {
    meta$cell_barcode <- as.character(meta$cell_barcode)
}

rownames(meta) <- meta$cell_barcode

colnames(scores)[1] <- "cell_barcode"
scores$cell_barcode <- as.character(scores$cell_barcode)
scores$Multinomial_Label <- as.character(scores$Multinomial_Label)

## =========================================================
## 2. MERGE
## =========================================================

common_cells <- Reduce(
    intersect,
    list(
        colnames(seurat_T),
        rownames(meta),
        scores$cell_barcode
    )
)

seurat_T <- subset(seurat_T, cells = common_cells)

meta_sub <- meta[common_cells, , drop = FALSE]
scores_sub <- scores[match(common_cells, scores$cell_barcode), , drop = FALSE]

for (cn in colnames(meta_sub)) {
    seurat_T@meta.data[[cn]] <- meta_sub[[cn]]
}

seurat_T$Multinomial_Label <- scores_sub$Multinomial_Label

cat("Cells after merge:", ncol(seurat_T), "\n")

## =========================================================
## 3. CREATE CAR GROUP
## =========================================================

seurat_T$CAR_group <- as.character(seurat_T$CAR_Variant)

seurat_T$CAR_group[
    seurat_T$condition == "Unstim_28z"
] <- "Unstim_CD28_CD3Z"

## =========================================================
## 4. FILTERING
## Keep:
## - stimulated Library CAR cells
## - Unstim CD28-CD3Z control cells
## =========================================================

keep_library <- (
    seurat_T$condition == "Library" &
        !is.na(seurat_T$CAR_Variant) &
        seurat_T$CAR_Variant != "" &
        seurat_T$CAR_Variant != "NA" &
        seurat_T$CAR_Variant != "D"
)

keep_unstim <- (
    seurat_T$condition == "Unstim_28z"
)

keep <- (
    (keep_library | keep_unstim) &
        !is.na(seurat_T$Multinomial_Label) &
        seurat_T$Multinomial_Label != ""
)

seurat_T <- subset(
    seurat_T,
    cells = colnames(seurat_T)[keep]
)

## clean labels
seurat_T$Multinomial_Label_clean <- gsub(
    "_",
    "-",
    seurat_T$Multinomial_Label
)

## remove technical labels
remove_labels <- c(
    "Doublet-RBC",
    "Doublet-Platelet",
    "Doublet-Myeloid",
    "Doublet-Bcell",
    "Doublet-Plasmablast",
    "Doublet-Fibroblast",
    "Poor-Quality"
)

seurat_T <- subset(
    seurat_T,
    cells = colnames(seurat_T)[
        !seurat_T$Multinomial_Label_clean %in% remove_labels
    ]
)

## =========================================================
## 5. FILTER CAR GROUPS BY CELL NUMBER
## Keep Unstim control regardless of threshold
## =========================================================

car_tab <- sort(
    table(seurat_T$CAR_group),
    decreasing = TRUE
)

valid_car <- names(car_tab)[
    car_tab >= 100 |
        names(car_tab) == "Unstim_CD28_CD3Z"
]

seurat_T <- subset(
    seurat_T,
    subset = CAR_group %in% valid_car
)

cat("Groups after filtering:\n")
print(sort(table(seurat_T$CAR_group), decreasing = TRUE))

car_counts_df <- data.frame(
    CAR_group = names(sort(table(seurat_T$CAR_group), decreasing = TRUE)),
    n_cells = as.integer(sort(table(seurat_T$CAR_group), decreasing = TRUE)),
    row.names = NULL
)

write.csv(
    car_counts_df,
    "Filtered_CAR_counts_with_Unstim_CD28_CD3Z.csv",
    row.names = FALSE
)

## =========================================================
## 6. CELL CYCLE SCORING
## =========================================================

seurat_T <- CellCycleScoring(
    object = seurat_T,
    s.features = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
)

cat("Cell cycle scoring done.\n")

## =========================================================
## 7. NKG2D LIGAND SCORE
## =========================================================

nkg2d_genes <- c(
    "MICA",
    "MICB",
    "ULBP1",
    "ULBP2",
    "ULBP3",
    "RAET1E",
    "RAET1G",
    "RAET1L"
)

nkg2d_found <- nkg2d_genes[
    nkg2d_genes %in% rownames(seurat_T)
]

cat("NKG2D ligand genes found:\n")
print(nkg2d_found)

if (length(nkg2d_found) == 0) {
    stop("No NKG2D ligand genes found in Seurat object.")
}

expr_mat <- GetAssayData(seurat_T, layer = "data")

if (length(nkg2d_found) == 1) {
    seurat_T$NKG2DL_score <- as.numeric(expr_mat[nkg2d_found, ])
} else {
    seurat_T$NKG2DL_score <- colMeans(
        expr_mat[nkg2d_found, , drop = FALSE],
        na.rm = TRUE
    )
}

## =========================================================
## 8. CELL-LEVEL TABLE
## =========================================================

plot_vars <- c(
    "S.Score",
    "G2M.Score",
    "NKG2DL_score"
)

cell_df <- FetchData(
    seurat_T,
    vars = plot_vars
)

cell_df$cell_barcode <- rownames(cell_df)
cell_df$CAR_group <- seurat_T$CAR_group
cell_df$condition <- seurat_T$condition
cell_df$Multinomial_Label <- seurat_T$Multinomial_Label_clean

write.csv(
    cell_df,
    "FilteredCAR_cell_table_CellCycle_NKG2DL_with_Unstim_CD28_CD3Z.csv",
    row.names = FALSE
)

## =========================================================
## 9. PIE COMPOSITION BY MULTINOMIAL LABEL
## =========================================================

pie_comp <- cell_df %>%
    group_by(CAR_group, Multinomial_Label) %>%
    summarise(
        n = n(),
        .groups = "drop"
    ) %>%
    group_by(CAR_group) %>%
    mutate(
        frac = n / sum(n)
    ) %>%
    ungroup()

write.csv(
    pie_comp,
    "FilteredCAR_pie_composition_with_Unstim_CD28_CD3Z.csv",
    row.names = FALSE
)

pie_wide <- pie_comp %>%
    tidyr::pivot_wider(
        names_from = Multinomial_Label,
        values_from = frac,
        values_fill = 0
    )

pie_cols <- setdiff(
    colnames(pie_wide),
    c("CAR_group", "n_cells", "n")
)

preferred_labels <- c(
    "CD4-CM",
    "CD4-EM",
    "CD4-Naive",
    "CD8-CM",
    "CD8-EM",
    "CD8-Naive",
    "CD8-TEMRA",
    "gdT",
    "MAIT",
    "Treg"
)

pie_cols <- c(
    preferred_labels[preferred_labels %in% pie_cols],
    setdiff(pie_cols, preferred_labels)
)

pie_wide <- pie_wide[, c("CAR_group", pie_cols), drop = FALSE]

## =========================================================
## 10. CAR COORDINATES
## mean per CAR/control group, then convert to Z-score
## =========================================================

car_xy <- cell_df %>%
    group_by(CAR_group) %>%
    summarise(
        n_cells = n(),
        S_mean = mean(S.Score, na.rm = TRUE),
        G2M_mean = mean(G2M.Score, na.rm = TRUE),
        NKG2DL_mean = mean(NKG2DL_score, na.rm = TRUE),
        .groups = "drop"
    )

car_xy$S_z <- as.numeric(scale(car_xy$S_mean))
car_xy$G2M_z <- as.numeric(scale(car_xy$G2M_mean))
car_xy$NKG2DL_z <- as.numeric(scale(car_xy$NKG2DL_mean))

plot_df <- car_xy %>%
    left_join(
        pie_wide,
        by = "CAR_group"
    )

plot_df$radius_fixed <- ifelse(
    plot_df$CAR_group == "Unstim_CD28_CD3Z",
    0.30,
    0.24
)

write.csv(
    plot_df,
    "FilteredCAR_PieScatter_coordinates_CellCycle_NKG2DL_Zscore_with_Unstim_CD28_CD3Z.csv",
    row.names = FALSE
)

## =========================================================
## 11. COLORS
## =========================================================

pie_fill <- c(
    "CD4-CM" = "#D55E5E",
    "CD4-EM" = "#C98A2E",
    "CD4-Naive" = "#B7A43A",
    "CD8-CM" = "#4DAF4A",
    "CD8-EM" = "#2FAE8F",
    "CD8-Naive" = "#3FA7C9",
    "CD8-TEMRA" = "#4F83C2",
    "gdT" = "#7B6FD0",
    "MAIT" = "#B96AC9",
    "Treg" = "#E07AAE"
)

missing_cols <- setdiff(
    pie_cols,
    names(pie_fill)
)

if (length(missing_cols) > 0) {
    extra_cols <- grDevices::rainbow(length(missing_cols))
    names(extra_cols) <- missing_cols
    pie_fill <- c(pie_fill, extra_cols)
}

pie_fill <- pie_fill[pie_cols]

## =========================================================
## 12. COMMON PLOT THEME
## =========================================================

plot_theme <- theme_classic(base_size = 18) +
    theme(
        axis.text = element_text(
            size = 18,
            color = "black"
        ),
        axis.title = element_text(
            size = 20,
            face = "bold",
            color = "black"
        ),
        plot.title = element_text(
            size = 22,
            face = "bold",
            hjust = 0.5,
            color = "black"
        ),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.2, "cm")
    )

label_offset <- 0.18

## =========================================================
## 13. PLOT 1 — S phase vs NKG2DL
## =========================================================

p_S <- ggplot() +
    geom_scatterpie(
        data = plot_df,
        aes(
            x = S_z,
            y = NKG2DL_z,
            r = radius_fixed
        ),
        cols = pie_cols,
        color = "black",
        size = 0.25,
        alpha = 0.95
    ) +
    geom_text(
        data = plot_df,
        aes(
            x = S_z,
            y = NKG2DL_z - label_offset,
            label = CAR_group
        ),
        size = 5.5,
        color = "black"
    ) +
    scale_fill_manual(values = pie_fill) +
    coord_equal() +
    plot_theme +
    labs(
        x = "S phase score (Z-score)",
        y = "NKG2DL score (Z-score)",
        title = "CAR T composition on S phase vs NKG2DL axes"
    )

ggsave(
    "FilteredCAR_PieScatter_S_vs_NKG2DL_with_labels_Zscore_with_Unstim_CD28_CD3Z.png",
    p_S,
    width = 14,
    height = 10,
    dpi = 300
)

p_S_nolabel <- ggplot() +
    geom_scatterpie(
        data = plot_df,
        aes(
            x = S_z,
            y = NKG2DL_z,
            r = radius_fixed
        ),
        cols = pie_cols,
        color = "black",
        size = 0.25,
        alpha = 0.95
    ) +
    scale_fill_manual(values = pie_fill) +
    coord_equal() +
    plot_theme +
    labs(
        x = "S phase score (Z-score)",
        y = "NKG2DL score (Z-score)",
        title = "CAR T composition on S phase vs NKG2DL axes"
    )

ggsave(
    "FilteredCAR_PieScatter_S_vs_NKG2DL_NO_LABELS_with_Unstim_CD28_CD3Z.png",
    p_S_nolabel,
    width = 14,
    height = 10,
    dpi = 300
)

## =========================================================
## 14. PLOT 2 — G2M vs NKG2DL
## =========================================================

p_G2M <- ggplot() +
    geom_scatterpie(
        data = plot_df,
        aes(
            x = G2M_z,
            y = NKG2DL_z,
            r = radius_fixed
        ),
        cols = pie_cols,
        color = "black",
        size = 0.25,
        alpha = 0.95
    ) +
    geom_text(
        data = plot_df,
        aes(
            x = G2M_z,
            y = NKG2DL_z - label_offset,
            label = CAR_group
        ),
        size = 5.5,
        color = "black"
    ) +
    scale_fill_manual(values = pie_fill) +
    coord_equal() +
    plot_theme +
    labs(
        x = "G2M score (Z-score)",
        y = "NKG2DL score (Z-score)",
        title = "CAR T composition on G2M vs NKG2DL axes"
    )

ggsave(
    "FilteredCAR_PieScatter_G2M_vs_NKG2DL_with_labels_Zscore_with_Unstim_CD28_CD3Z.png",
    p_G2M,
    width = 14,
    height = 10,
    dpi = 300
)

p_G2M_nolabel <- ggplot() +
    geom_scatterpie(
        data = plot_df,
        aes(
            x = G2M_z,
            y = NKG2DL_z,
            r = radius_fixed
        ),
        cols = pie_cols,
        color = "black",
        size = 0.25,
        alpha = 0.95
    ) +
    scale_fill_manual(values = pie_fill) +
    coord_equal() +
    plot_theme +
    labs(
        x = "G2M score (Z-score)",
        y = "NKG2DL score (Z-score)",
        title = "CAR T composition on G2M vs NKG2DL axes"
    )

ggsave(
    "FilteredCAR_PieScatter_G2M_vs_NKG2DL_NO_LABELS_with_Unstim_CD28_CD3Z.png",
    p_G2M_nolabel,
    width = 14,
    height = 10,
    dpi = 300
)

cat("\nDone.\n")
cat("Created files with Unstim_CD28_CD3Z control.\n")
