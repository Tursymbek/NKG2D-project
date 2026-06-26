# -----------------------------------------------------------------------------
# GSE214231 Seurat preprocessing and starCAT export
#
# Original archived script: Castellanos-Rueda et al. [76]/1 - Data preparation TCAT.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

## =========================================================
## starCAT (data preparation and launching starCAT)
## =========================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(readr)
    library(R.utils)
    library(ggplot2)
    library(dplyr)
})

## =========================================================
## 1. DIRECTORIES 
## =========================================================
base_path <- Sys.getenv("GSE214231_RAW_DIR", unset = "data/external/GSE214231_RAW")
out_dir   <- file.path(base_path, "starCAT_run")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. FILE DESCRIPTION
## =========================================================
sample_table <- data.frame(
    sample_name = c(
        "D1_2DomLib",
        "D1_2DomLib_WT",
        "D2_2DomLib_WT",
        "D2_Unstim_28z",
        "D3_2DomLib_1",
        "D3_2DomLib_2",
        "D3_2DomLib_WT",
        "D3_hT"
    ),
    folder_name = c(
        "D1_2DomLib",
        "D1_2DomLib_WT",
        "D2_2DomLib_WT",
        "D2_Unstim_28z",
        "D3_2DomLib_1",
        "D3_2DomLib_2",
        "D3_2DomLib_WT",
        "D3_hT"
    ),
    donor_id = c("D1","D1","D2","D2","D3","D3","D3","D3"),
    condition = c("Library","WT","WT","Unstim_28z","Library","Library","WT","hT"),
    stringsAsFactors = FALSE
)

barcode_map <- c(
    "D1_2DomLib"    = "GSM6601767_D1_2DomLib_barcodes_assigned.csv",
    "D1_2DomLib_WT" = "GSM6601768_D1_2DomLib_WT_barcodes_assigned.csv",
    "D2_2DomLib_WT" = "GSM6601769_D2_2DomLib_WT_barcodes_assigned.csv",
    "D3_2DomLib_1"  = "GSM6601770_D3_2DomLib_1_barcodes_assigned.csv",
    "D3_2DomLib_2"  = "GSM6601771_D3_2DomLib_2_barcodes_assigned.csv",
    "D3_2DomLib_WT" = "GSM6601772_D3_2DomLib_WT_barcodes_assigned.csv"
)

read_barcode_csv <- function(csv_path) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
    
    out <- data.frame(
        barcode = trimws(as.character(df[["barcodes_assigned"]])),
        CAR_Variant = as.character(df[["CAR"]]),
        stringsAsFactors = FALSE
    )
    
    out <- out[!is.na(out$barcode) & out$barcode != "", , drop = FALSE]
    out$barcode <- ifelse(grepl("-1$", out$barcode), out$barcode, paste0(out$barcode, "-1"))
    out
}


## NOTE:
## This preparation script intentionally keeps all samples:
## - Library: CAR library after tumor co-culture
## - WT / Unstim_28z / hT: controls
## starCAT annotation is run on the combined T-cell object.
## Downstream CAR-ranking scripts should filter to condition == "Library".

## =========================================================
## 3. READING ALL 10X FOLDERS
## =========================================================
seu_list <- list()

for (i in seq_len(nrow(sample_table))) {
    s <- sample_table$sample_name[i]
    folder <- file.path(base_path, sample_table$folder_name[i])
    
    cat("Reading:", s, "\n")
    
    counts <- Read10X(folder)
    seu <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
    seu <- RenameCells(seu, add.cell.id = s)
    
    seu$sample_name <- s
    seu$donor_id <- sample_table$donor_id[i]
    seu$condition <- sample_table$condition[i]
    seu$CAR_Variant <- NA_character_
    
    if (s %in% names(barcode_map)) {
        ann <- read_barcode_csv(file.path(base_path, barcode_map[[s]]))
        ann$merged_barcode <- paste0(s, "_", ann$barcode)
        idx <- match(colnames(seu), ann$merged_barcode)
        seu$CAR_Variant <- ann$CAR_Variant[idx]
    }
    
    seu_list[[s]] <- seu
}

## =========================================================
## 4. MERGE
## =========================================================
seurat_all <- Reduce(function(x, y) merge(x, y = y), seu_list)
DefaultAssay(seurat_all) <- "RNA"

## =========================================================
## 5. QC
## =========================================================
seurat_all[["percent.mt"]] <- PercentageFeatureSet(seurat_all, pattern = "^MT-")

seurat_all <- subset(
    seurat_all,
    subset = nCount_RNA >= 500 &
        nFeature_RNA >= 200 &
        nFeature_RNA <= 7000 &
        percent.mt <= 15
)

## =========================================================
## 6. NORMALIZATION FOR T CELL FILTERING
## =========================================================
seurat_all <- NormalizeData(seurat_all, verbose = FALSE)
seurat_all <- FindVariableFeatures(seurat_all, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seurat_all <- ScaleData(seurat_all, verbose = FALSE)
seurat_all <- RunPCA(seurat_all, npcs = 30, verbose = FALSE)
seurat_all <- RunUMAP(seurat_all, dims = 1:20, verbose = FALSE)
seurat_all <- FindNeighbors(seurat_all, dims = 1:20, verbose = FALSE)
seurat_all <- FindClusters(seurat_all, resolution = 0.4, verbose = FALSE)

## =========================================================
## 7.T CELLS ALLOCATION
## =========================================================
get_gene_or_zero <- function(obj, gene) {
    if (gene %in% rownames(obj)) {
        FetchData(obj, vars = gene)[,1]
    } else {
        rep(0, ncol(obj))
    }
}

seurat_all$PTPRC_expr <- get_gene_or_zero(seurat_all, "PTPRC")
seurat_all$CD3D_expr  <- get_gene_or_zero(seurat_all, "CD3D")
seurat_all$CD3E_expr  <- get_gene_or_zero(seurat_all, "CD3E")
seurat_all$EPCAM_expr <- get_gene_or_zero(seurat_all, "EPCAM")
seurat_all$ERBB2_expr <- get_gene_or_zero(seurat_all, "ERBB2")
seurat_all$KRT8_expr  <- get_gene_or_zero(seurat_all, "KRT8")
seurat_all$KRT18_expr <- get_gene_or_zero(seurat_all, "KRT18")

tcell_features <- intersect(c("PTPRC","CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB","NKG7","CD4","CD8A"), rownames(seurat_all))
tumor_features <- intersect(c("EPCAM","KRT8","KRT18","KRT19","ERBB2"), rownames(seurat_all))

if (length(tcell_features) > 0) {
    seurat_all <- AddModuleScore(seurat_all, features = list(tcell_features), name = "TcellScore")
} else {
    seurat_all$TcellScore1 <- 0
}

if (length(tumor_features) > 0) {
    seurat_all <- AddModuleScore(seurat_all, features = list(tumor_features), name = "TumorScore")
} else {
    seurat_all$TumorScore1 <- 0
}

seurat_T <- subset(
    seurat_all,
    subset =
        (
            !is.na(CAR_Variant) |
                (
                    (PTPRC_expr > 0 | CD3D_expr > 0 | CD3E_expr > 0 | TcellScore1 > 0) &
                        TumorScore1 <= TcellScore1
                )
        ) &
        EPCAM_expr == 0 &
        ERBB2_expr == 0 &
        KRT8_expr == 0 &
        KRT18_expr == 0
)

## =========================================================
## 8. UMAP RECALCULATION TO T CELLS
## =========================================================
seurat_T <- NormalizeData(seurat_T, verbose = FALSE)
seurat_T <- FindVariableFeatures(seurat_T, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seurat_T <- ScaleData(seurat_T, verbose = FALSE)
seurat_T <- RunPCA(seurat_T, npcs = 30, verbose = FALSE)
seurat_T <- RunUMAP(seurat_T, dims = 1:20, verbose = FALSE)
seurat_T <- FindNeighbors(seurat_T, dims = 1:20, verbose = FALSE)
seurat_T <- FindClusters(seurat_T, resolution = 0.4, verbose = FALSE)

saveRDS(seurat_T, file.path(out_dir, "seurat_T.rds"))

## =========================================================
## 9. EXPORT TO STARCAT FORMAT
## =========================================================
DefaultAssay(seurat_T) <- "RNA"

## важно для Seurat v5 после merge
seurat_T <- JoinLayers(seurat_T, assay = "RNA")

## raw counts
counts <- LayerData(seurat_T, assay = "RNA", layer = "counts")

keep_genes <- Matrix::rowSums(counts > 0) >= 10
counts <- counts[keep_genes, , drop = FALSE]

writeMM(counts, file.path(out_dir, "matrix.mtx"))
gzip(file.path(out_dir, "matrix.mtx"), overwrite = TRUE)

write_delim(
  data.frame(barcodes = colnames(counts)),
  file.path(out_dir, "barcodes.tsv"),
  delim = "\t",
  col_names = FALSE
)
gzip(file.path(out_dir, "barcodes.tsv"), overwrite = TRUE)

features <- data.frame(
  gene_id = rownames(counts),
  gene_name = rownames(counts),
  type = "Gene Expression"
)
write_delim(
  features,
  file.path(out_dir, "features.tsv"),
  delim = "\t",
  col_names = FALSE
)
gzip(file.path(out_dir, "features.tsv"), overwrite = TRUE)

## =========================================================
## 10. METADATA SAVING
## =========================================================
metadata_cells <- seurat_T@meta.data
metadata_cells$cell_barcode <- colnames(seurat_T)
write.csv(metadata_cells, file.path(out_dir, "metadata_cells.csv"), row.names = FALSE)
write.csv(sample_table, file.path(out_dir, "metadata_samples.csv"), row.names = FALSE)

## =========================================================
## 11. STARCAT LAUNCHING
## =========================================================
output_name <- "GSE214231_Tcells"

starcat_exe <- "C:/Users/Lenovo/AppData/Local/Packages/PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0/LocalCache/local-packages/Python310/Scripts/starcat.exe"

cmd <- paste0(
    '"', starcat_exe, '"',
    ' --reference "TCAT.V1"',
    ' --counts "', file.path(out_dir, "matrix.mtx.gz"), '"',
    ' --output-dir "', out_dir, '"',
    ' --name "', output_name, '"'
)

cat(cmd, "\n")
system(cmd)

## =========================================================
## 12. READING STARCAT RESULTS
## =========================================================
usage_file  <- file.path(out_dir, paste0(output_name, ".rf_usage_normalized.txt"))
scores_file <- file.path(out_dir, paste0(output_name, ".scores.txt"))

if (file.exists(usage_file)) {
    usage <- read.table(usage_file, header = TRUE, sep = "\t", check.names = FALSE)
    usage_barcodes <- usage[,1]
    usage_mat <- usage[,-1, drop = FALSE]
    rownames(usage_mat) <- usage_barcodes
    
    common_cells <- intersect(colnames(seurat_T), rownames(usage_mat))
    seurat_T <- subset(seurat_T, cells = common_cells)
    usage_mat <- usage_mat[common_cells, , drop = FALSE]
    
    for (nm in colnames(usage_mat)) {
        seurat_T[[paste0("TCAT_", nm)]] <- usage_mat[, nm]
    }
    
    seurat_T$TCAT_top_program <- colnames(usage_mat)[max.col(usage_mat, ties.method = "first")]
    
    saveRDS(seurat_T, file.path(out_dir, "seurat_T_with_starCAT.rds"))
    write.csv(seurat_T@meta.data, file.path(out_dir, "metadata_cells_with_starCAT.csv"))
    
    cat("starCAT results loaded successfully.\n")
} else {
    cat("starCAT output not found yet. Check installation or run the printed command manually.\n")
}
