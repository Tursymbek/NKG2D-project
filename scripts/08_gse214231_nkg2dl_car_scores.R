# -----------------------------------------------------------------------------
# GSE214231 CAR-level NKG2DL score analysis
#
# Original archived script: Castellanos-Rueda et al. [76]/B - NKG2D ligand SUM score by CAR Library.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

## =========================================================
## NKG2D ligand SUM score by CAR
## Includes Unstim CD28-CD3Z control
## CD223-CD79B red, Unstim CD28-CD3Z blue
## =========================================================

suppressPackageStartupMessages({
    library(Seurat)    
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(ggtext)
})

## =========================================================
## 1. READ DATA
## =========================================================

meta <- read.csv(
    "metadata_cells_with_starCAT.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (colnames(meta)[1] == "") {
    colnames(meta)[1] <- "cell_barcode"
}

## =========================================================
## 2. DEFINE NKG2D LIGANDS
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

nice_names <- c(
    "MICA"   = "MICA",
    "MICB"   = "MICB",
    "ULBP1"  = "ULBP1",
    "ULBP2"  = "ULBP2",
    "ULBP3"  = "ULBP3",
    "RAET1E" = "ULBP4",
    "RAET1G" = "ULBP5",
    "RAET1L" = "ULBP6"
)

## =========================================================
## 3. LOAD SEURAT OBJECT
## =========================================================

seu <- readRDS("seurat_T_with_starCAT.rds")

DefaultAssay(seu) <- "RNA"

expr <- FetchData(
    seu,
    vars = nkg2d_genes
)

expr$cell_barcode <- rownames(expr)

## =========================================================
## 4. MERGE WITH METADATA
## =========================================================

meta$cell_barcode <- trimws(as.character(meta$cell_barcode))
expr$cell_barcode <- trimws(as.character(expr$cell_barcode))

dat <- left_join(
    meta,
    expr,
    by = "cell_barcode"
)

## =========================================================
## 5. FILTER: Library CAR cells + Unstim CD28-CD3Z control
## =========================================================

dat <- dat %>%
    filter(
        (
            condition == "Library" &
                !is.na(CAR_Variant) &
                CAR_Variant != "" &
                CAR_Variant != "NA" &
                CAR_Variant != "D"
        ) |
            condition == "Unstim_28z"
    )

cat("Cells used before CAR filtering:", nrow(dat), "\n")
print(table(dat$condition, useNA = "ifany"))
print(table(dat$sample_name, useNA = "ifany"))

car_counts <- dat %>%
    filter(condition == "Library") %>%
    count(CAR_Variant, name = "n_cells")

valid_cars <- car_counts %>%
    filter(n_cells >= 100) %>%
    pull(CAR_Variant)

dat <- dat %>%
    filter(
        (
            condition == "Library" &
                CAR_Variant %in% valid_cars
        ) |
            condition == "Unstim_28z"
    )

dat$CAR_plot <- ifelse(
    dat$condition == "Unstim_28z",
    "Unstim CD28-CD3Z",
    as.character(dat$CAR_Variant)
)

cat("Cells used after CAR filtering:", nrow(dat), "\n")
print(table(dat$CAR_plot, useNA = "ifany"))

## =========================================================
## 6. CALCULATE TOTAL NKG2DL SCORE
## SUM instead of MEAN
## =========================================================

present_genes <- intersect(nkg2d_genes, colnames(dat))

cat("NKG2D ligand columns found in dat:\n")
print(present_genes)

if (length(present_genes) == 0) {
    stop("No NKG2D ligand columns found in dat. Check colnames(dat).")
}

dat$NKG2DL_sum <- rowSums(
    dat[, present_genes, drop = FALSE],
    na.rm = TRUE
)

## =========================================================
## 7. SUM SCORE PER CAR / CONTROL
## =========================================================

plot_df <- dat %>%
    group_by(CAR_plot) %>%
    summarise(
        NKG2DL_sum = mean(NKG2DL_sum, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
    ) %>%
    arrange(desc(NKG2DL_sum))

write.csv(
    plot_df,
    "Library_plus_Unstim28z_CAR_NKG2DL_sum_scores.csv",
    row.names = FALSE
)

## =========================================================
## 8. COLORS
## =========================================================

cellpress_blue <- "#4A90C2"
highlight_red  <- "#E64B35"
other_gray     <- "#BFBFBF"

plot_df$bar_color <- dplyr::case_when(
    plot_df$CAR_plot == "CD223-CD79B" ~ highlight_red,
    plot_df$CAR_plot == "Unstim CD28-CD3Z" ~ cellpress_blue,
    TRUE ~ other_gray
)

plot_df$CAR_label <- dplyr::case_when(
    plot_df$CAR_plot == "CD223-CD79B" ~
        paste0("<span style='color:", highlight_red, ";'>CD223-CD79B</span>"),
    
    plot_df$CAR_plot == "Unstim CD28-CD3Z" ~
        paste0("<span style='color:", cellpress_blue, ";'>Unstim CD28-CD3Z</span>"),
    
    TRUE ~ plot_df$CAR_plot
)

plot_df$CAR_label <- factor(
    plot_df$CAR_label,
    levels = rev(plot_df$CAR_label)
)

## =========================================================
## 9. PLOT
## =========================================================

p <- ggplot(
    plot_df,
    aes(
        x = NKG2DL_sum,
        y = CAR_label
    )
) +
    geom_col(
        aes(fill = bar_color),
        width = 0.85
    ) +
    scale_fill_identity() +
    theme_classic(base_size = 18) +
    theme(
        axis.text.x = element_text(
            size = 18,
            color = "black"
        ),
        axis.text.y = ggtext::element_markdown(
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
            hjust = 0.5
        )
    ) +
    labs(
        title = "Total NKG2D ligand score by CAR",
        x = "Total NKG2D ligand score",
        y = "CAR"
    )

ggsave(
    "Library_plus_Unstim28z_CAR_total_NKG2DL_score.png",
    p,
    width = 10,
    height = 7,
    dpi = 300
)

ggsave(
    "Library_plus_Unstim28z_CAR_total_NKG2DL_score.pdf",
    p,
    width = 10,
    height = 7
)

p
