# =========================================================
# FILTER TCGA SAMPLES BASED ON EXCEL LIST (KEEP GTEX)
# =========================================================

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(dplyr)
})

# =========================================================
# 1. FILE PATHS
# =========================================================
expr_file  <- "TcgaTargetGtex_rsem_gene_tpm"
excel_file <- "41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"

# =========================================================
# 2. LOAD DATA
# =========================================================
cat("Loading expression...\n")
expr <- fread(expr_file, data.table = FALSE)

cat("Loading excel...\n")
meta <- read_excel(excel_file)

# =========================================================
# 3. CHECK COLUMN WITH SAMPLE IDS
# =========================================================
cat("Excel columns:\n")
print(colnames(meta))

# УКАЖИ правильную колонку, если имя отличается
sample_col <- "Sample ID"

if (!sample_col %in% colnames(meta)) {
    stop(paste0("Column '", sample_col, "' not found in Excel. Check colnames(meta)."))
}

sample_ids <- as.character(meta[[sample_col]])
sample_ids <- trimws(sample_ids)
sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]

cat("Excel sample count:", length(sample_ids), "\n")
cat("Example Excel IDs:\n")
print(head(sample_ids, 10))

# =========================================================
# 4. EXTRACT SAMPLE NAMES FROM EXPRESSION MATRIX
# =========================================================
expr_samples <- colnames(expr)[-1]

cat("Expression sample count:", length(expr_samples), "\n")
cat("Example expression IDs:\n")
print(head(expr_samples, 10))

# =========================================================
# 5. DEFINE TCGA vs GTEX
# =========================================================
is_tcga <- grepl("^TCGA", expr_samples)
is_gtex <- grepl("^GTEX", expr_samples)

tcga_samples <- expr_samples[is_tcga]
gtex_samples <- expr_samples[is_gtex]

cat("TCGA samples:", length(tcga_samples), "\n")
cat("GTEx samples:", length(gtex_samples), "\n")

# =========================================================
# 6. NORMALIZE IDS AND FILTER TCGA
# =========================================================
# Приводим разделители к одному виду
sample_ids  <- gsub("\\.", "-", sample_ids)
sample_ids  <- gsub("_", "-", sample_ids)

tcga_samples_clean <- gsub("\\.", "-", tcga_samples)
tcga_samples_clean <- gsub("_", "-", tcga_samples_clean)

# Сравниваем по patient barcode: TCGA-XX-XXXX (12 символов)
short_tcga <- substr(tcga_samples_clean, 1, 12)
short_meta <- substr(sample_ids,        1, 12)

cat("Unique TCGA patient IDs in expression:", length(unique(short_tcga)), "\n")
cat("Unique TCGA patient IDs in Excel:", length(unique(short_meta)), "\n")
cat("Matching TCGA patient IDs:", sum(unique(short_tcga) %in% unique(short_meta)), "\n")

# Оставляем TCGA samples, чьи patient IDs есть в Excel
tcga_keep <- tcga_samples[short_tcga %in% short_meta]

cat("TCGA kept:", length(tcga_keep), "\n")
cat("Example kept TCGA IDs:\n")
print(head(tcga_keep, 10))

if (length(tcga_keep) == 0) {
    stop("No TCGA samples matched. Check sample_col and barcode format in Excel.")
}

# =========================================================
# 7. FINAL SAMPLE LIST
# =========================================================
final_samples <- c(tcga_keep, gtex_samples)

cat("Final sample count:", length(final_samples), "\n")

# =========================================================
# 8. SUBSET EXPRESSION MATRIX
# =========================================================
expr_filtered <- expr[, c(1, which(colnames(expr) %in% final_samples)), drop = FALSE]

cat("Filtered expression dim:", dim(expr_filtered), "\n")

# =========================================================
# 9. SAVE
# =========================================================
fwrite(expr_filtered, "TcgaTargetGtex_filtered.tsv", sep = "\t", quote = FALSE)

cat("DONE: saved as TcgaTargetGtex_filtered.tsv\n")
