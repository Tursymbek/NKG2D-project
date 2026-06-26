# -----------------------------------------------------------------------------
# GSE214231 clean correlation analysis with NKG2D ligands
#
# Original archived script: Castellanos-Rueda et al. [76]/06_Correlation_with_ligands_clean.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

## =========================================================
## Correlation of TCAT T cell program and NKG2DL
## =========================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tibble)
})



## =========================================================
## 1. LOAD SEURAT
## =========================================================
seurat_T <- readRDS("seurat_T_with_starCAT.rds")

meta <- seurat_T@meta.data
meta$cell_barcode <- rownames(meta)

## =========================================================
## FILTER: only stimulated CAR library cells
## =========================================================
## Correlation is calculated only within CAR library cells after tumor co-culture.
## WT, Unstim_28z and hT controls are excluded from CAR-variant interpretation.

meta <- meta %>%
    filter(
        condition == "Library",
        !is.na(CAR_Variant),
        CAR_Variant != "",
        CAR_Variant != "NA",
        CAR_Variant != "D"
    )

cat("Library CAR cells used:", nrow(meta), "\n")
print(table(meta$sample_name, useNA = "ifany"))

## =========================================================
## 2. DEFINE T-CELL PROGRAMS TO KEEP
## =========================================================
remove_exact <- c(
    "Doublet-RBC",
    "Doublet-Platelet",
    "Doublet-Myeloid",
    "Doublet-Bcell",
    "Doublet-Plasmablast",
    "Doublet-Fibroblast",
    "Poor-Quality",
    "CD172a/MERTK",
    "CD172a.MERTK"
)

keep_programs <- c(
    "CD4-CM", "CD4-Naive",
    "CD8-EM", "CD8-Naive", "CD8-Trm", "TEMRA",
    "Treg", "MAIT", "gdT",
    "Th1-Like", "Th1-1", "Th2", "Th2-Activated", "Th2-Resting",
    "Th17-Activated", "Th17-Resting", "Th22",
    "Tfh-1", "Tfh-2", "Tph",
    "Exhaustion", "Cytotoxic",
    "ICOS/CD38", "CTLA4/CD38", "OX40/EBI3", "TIMD4/TIM3",
    "SOX4/TOX2", "CD40LG/TXNIP", "IL10/IL19",
    "ISG", "HLA", "IEG", "IEG2", "IEG3",
    "CellCycle-S", "CellCycle-G2M", "CellCycle-Late-S"
)

## =========================================================
## 3. FIND AVAILABLE TCAT COLUMNS
## =========================================================
tcat_cols <- grep("^TCAT_", colnames(meta), value = TRUE)

## clean names
program_map <- data.frame(
    Program = tcat_cols,
    Program_clean = sub("^TCAT_", "", tcat_cols),
    stringsAsFactors = FALSE
)

program_map <- program_map %>%
    filter(!Program_clean %in% remove_exact) %>%
    filter(!grepl("^Doublet", Program_clean)) %>%
    filter(Program_clean %in% keep_programs)

print(program_map)

## =========================================================
## 4. LIGANDS + TOTAL NKG2DL SCORE
## =========================================================
ligands <- c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L")
ligands_found <- ligands[ligands %in% rownames(seurat_T)]

cat("Ligands found:\n")
print(ligands_found)

if (length(ligands_found) == 0) {
    stop("No NKG2D ligands found in Seurat object.")
}

expr_df <- FetchData(seurat_T, vars = ligands_found, cells = meta$cell_barcode)
expr_df$cell_barcode <- rownames(expr_df)

## true total ligand score = SUM across available genes
expr_df$NKG2DL_sum <- rowSums(expr_df[, ligands_found, drop = FALSE], na.rm = TRUE)

## =========================================================
## 5. MERGE EXPRESSION + TCAT PROGRAM USAGE
## =========================================================
dat <- meta %>%
    select(cell_barcode, all_of(program_map$Program)) %>%
    inner_join(expr_df, by = "cell_barcode")

cat("Merged cells:", nrow(dat), "\n")

## =========================================================
## 6. COMPUTE CORRELATIONS
## =========================================================
genes_to_test <- c(ligands_found, "NKG2DL_sum")

corr_list <- list()

for (i in seq_len(nrow(program_map))) {
    prog_col   <- program_map$Program[i]
    prog_clean <- program_map$Program_clean[i]
    
    for (g in genes_to_test) {
        tmp <- dat[, c(prog_col, g)]
        tmp <- tmp[complete.cases(tmp), , drop = FALSE]
        
        if (nrow(tmp) < 30) next
        if (length(unique(tmp[[prog_col]])) < 2) next
        if (length(unique(tmp[[g]])) < 2) next
        
        ct <- suppressWarnings(cor.test(
            tmp[[prog_col]],
            tmp[[g]],
            method = "spearman",
            exact = FALSE
        ))
        
        corr_list[[length(corr_list) + 1]] <- data.frame(
            Program = prog_col,
            Program_clean = prog_clean,
            Gene = g,
            rho = unname(ct$estimate),
            p_value = ct$p.value,
            n_cells = nrow(tmp),
            stringsAsFactors = FALSE
        )
    }
}

corr <- bind_rows(corr_list)

## FDR
corr$FDR <- p.adjust(corr$p_value, method = "BH")

write.csv(corr, "Library_only_TCAT_vs_NKG2DL_correlation_table_with_SUM.csv", row.names = FALSE)

## =========================================================
## 7. ORDER PROGRAMS AND RENAME LIGANDS
## =========================================================

program_order <- c(
    "CD4-Naive", "CD4-CM",
    "CD8-Naive", "CD8-EM", "CD8-Trm", "TEMRA",
    "Treg", "MAIT", "gdT",
    "Th1-Like", "Th1-1",
    "Th2", "Th2-Resting", "Th2-Activated",
    "Th17-Resting", "Th17-Activated", "Th22",
    "Tfh-1", "Tfh-2", "Tph",
    "Cytotoxic", "Exhaustion",
    "ICOS/CD38", "CTLA4/CD38", "OX40/EBI3", "TIMD4/TIM3",
    "SOX4/TOX2", "CD40LG/TXNIP", "IL10/IL19",
    "IEG", "IEG2", "IEG3",
    "ISG", "HLA",
    "CellCycle-S", "CellCycle-Late-S", "CellCycle-G2M"
)

## Rename ligands
gene_labels <- c(
    "MICA" = "MICA",
    "MICB" = "MICB",
    "ULBP1" = "ULBP1",
    "ULBP2" = "ULBP2",
    "ULBP3" = "ULBP3",
    "RAET1E" = "ULBP4",
    "RAET1G" = "ULBP5",
    "RAET1L" = "ULBP6",
    "NKG2DL_sum" = "Total NKG2DL"
)

gene_order <- names(gene_labels)

corr$Program_clean <- factor(corr$Program_clean,
                             levels = rev(program_order))

corr$Gene <- factor(corr$Gene,
                    levels = gene_order)

corr_plot <- corr %>%
    filter(!is.na(Program_clean),
           !is.na(Gene))

## =========================================================
## 8. HEATMAP
## =========================================================

p1 <- ggplot(corr_plot,
             aes(x = Gene,
                 y = Program_clean,
                 fill = rho)) +
    
    geom_tile(color = "grey85",
              linewidth = 0.4) +
    
    scale_x_discrete(labels = gene_labels) +
    
    scale_fill_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        na.value = "grey90"
    ) +
    
    theme_classic(base_size = 18) +
    
    theme(
        axis.text.x = element_text(
            size = 18,
            angle = 45,
            hjust = 1,
            face = "bold"
        ),
        
        axis.text.y = element_text(
            size = 18
        ),
        
        axis.title.x = element_text(
            size = 22,
            face = "bold"
        ),
        
        axis.title.y = element_text(
            size = 22,
            face = "bold"
        ),
        
        plot.title = element_text(
            size = 24,
            face = "bold",
            hjust = 0.5
        ),
        
        legend.title = element_text(
            size = 18,
            face = "bold"
        ),
        
        legend.text = element_text(
            size = 16
        )
    ) +
    
    labs(
        title = "Correlation between TCAT programs and NKG2D ligands",
        x = "NKG2D ligand",
        y = "T-cell TCAT program",
        fill = "Spearman rho"
    )

ggsave(
    "Library_only_TCAT_vs_NKG2DL_correlation_heatmap_with_SUM.png",
    p1,
    width = 12,
    height = 14,
    dpi = 300
)

## =========================================================
## 9. DOTPLOT
## =========================================================

corr_plot$logFDR <- -log10(corr_plot$FDR)
corr_plot$logFDR[!is.finite(corr_plot$logFDR)] <- NA

p2 <- ggplot(corr_plot,
             aes(x = Gene,
                 y = Program_clean,
                 size = logFDR,
                 color = rho)) +
    
    geom_point() +
    
    scale_x_discrete(labels = gene_labels) +
    
    scale_color_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        na.value = "grey80"
    ) +
    
    theme_classic(base_size = 18) +
    
    theme(
        axis.text.x = element_text(
            size = 18,
            angle = 45,
            hjust = 1,
            face = "bold"
        ),
        
        axis.text.y = element_text(
            size = 18
        ),
        
        axis.title.x = element_text(
            size = 22,
            face = "bold"
        ),
        
        axis.title.y = element_text(
            size = 22,
            face = "bold"
        ),
        
        plot.title = element_text(
            size = 24,
            face = "bold",
            hjust = 0.5
        ),
        
        legend.title = element_text(
            size = 18,
            face = "bold"
        ),
        
        legend.text = element_text(
            size = 16
        )
    ) +
    
    labs(
        title = "TCAT programs vs NKG2D ligands",
        x = "NKG2D ligand",
        y = "T-cell TCAT program",
        color = "Spearman rho",
        size = "-log10(FDR)"
    )

ggsave(
    "Library_only_TCAT_vs_NKG2DL_correlation_dotplot_with_SUM.png",
    p2,
    width = 12,
    height = 14,
    dpi = 300
)

## =========================================================
## 10. TOP HITS
## =========================================================
top_hits <- corr_plot %>%
    filter(!is.na(rho), !is.na(FDR)) %>%
    arrange(FDR, desc(abs(rho)))

write.csv(top_hits, "Library_only_TCAT_vs_NKG2DL_top_hits_with_SUM.csv", row.names = FALSE)

cat("Done.\n")
cat("Created files:\n")
cat("- Library_only_TCAT_vs_NKG2DL_correlation_table_with_SUM.csv\n")
cat("- Library_only_TCAT_vs_NKG2DL_correlation_heatmap_with_SUM.png\n")
cat("- Library_only_TCAT_vs_NKG2DL_correlation_dotplot_with_SUM.png\n")
cat("- Library_only_TCAT_vs_NKG2DL_top_hits_with_SUM.csv\n")
