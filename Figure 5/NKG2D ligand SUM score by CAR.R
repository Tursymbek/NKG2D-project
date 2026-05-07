## =========================================================
## NKG2D ligand SUM score by CAR
## Cell Press-style blue
## =========================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
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

## rename for plotting
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
## 5. FILTER CARS
## =========================================================

car_counts <- dat %>%
    filter(
        !is.na(CAR_Variant),
        CAR_Variant != "",
        CAR_Variant != "NA",
        CAR_Variant != "D"
    ) %>%
    count(CAR_Variant, name = "n_cells")

valid_cars <- car_counts %>%
    filter(n_cells >= 100) %>%
    pull(CAR_Variant)

dat <- dat %>%
    filter(CAR_Variant %in% valid_cars)

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
## 7. SUM SCORE PER CAR
## =========================================================

plot_df <- dat %>%
    group_by(CAR_Variant) %>%
    summarise(
        NKG2DL_sum = mean(NKG2DL_sum, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
    ) %>%
    arrange(desc(NKG2DL_sum))

plot_df$CAR_Variant <- factor(
    plot_df$CAR_Variant,
    levels = rev(plot_df$CAR_Variant)
)

write.csv(
    plot_df,
    "CAR_NKG2DL_sum_scores.csv",
    row.names = FALSE
)

## =========================================================
## 8. CELL PRESS-LIKE BLUE
## =========================================================

cellpress_blue <- "#4A90C2"

## =========================================================
## 9. PLOT
## =========================================================

p <- ggplot(
    plot_df,
    aes(
        x = NKG2DL_sum,
        y = CAR_Variant
    )
) +
    geom_col(
        fill = cellpress_blue,
        width = 0.85
    ) +
    theme_classic(base_size = 18) +
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
            hjust = 0.5
        )
    ) +
    labs(
        title = "Total NKG2D ligand score by CAR",
        x = "Total NKG2D ligand score",
        y = "CAR"
    )

ggsave(
    "CAR_total_NKG2DL_score.png",
    p,
    width = 10,
    height = 7,
    dpi = 300
)

p
