library(data.table)

# =========================================================
# 1. FILES
# =========================================================
expr_file  <- "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
model_file <- "Model.csv"

expr  <- fread(expr_file)
model <- fread(model_file)

# =========================================================
# 2. TARGET CELL LINES
# Use StrippedCellLineName from Model.csv
# =========================================================
target_order <- data.table(
    CellLineName = c(
        "HCT 116",
        "A549",
        "PC-3",
        "K-562",
        "SK-N-SH",
        "NCI-H929",
        "Hep G2",
        "MCF7",
        "PANC-1"
    ),
    StrippedCellLineName = c(
        "HCT116",
        "A549",
        "PC3",
        "K562",
        "SKNSH",
        "NCIH929",
        "HEPG2",
        "MCF7",
        "PANC1"
    ),
    order_id = 1:9
)

# =========================================================
# 3. FIND MODEL IDs USING Model.csv
# =========================================================
tmp <- merge(
    target_order,
    model[, .(ModelID, CellLineName, StrippedCellLineName, CCLEName)],
    by = "StrippedCellLineName",
    all.x = TRUE,
    suffixes = c("_target", "_DepMap")
)

setorder(tmp, order_id)

cat("\nMatched cell lines:\n")
print(tmp)

cat("\nNot found in Model.csv:\n")
print(tmp[is.na(ModelID)])

# stop if any line was not found
if (any(is.na(tmp$ModelID))) {
    stop("Some cell lines were not found in Model.csv")
}

# =========================================================
# 4. FILTER EXPRESSION BY ModelID
# In this DepMap file:
# rows = cell lines
# columns = genes
# =========================================================
expr_sub <- expr[ModelID %in% tmp$ModelID]

cat("\nExpression rows found:\n")
print(expr_sub[, .(ModelID, SequencingID, IsDefaultEntryForModel)])

cat("\nMissing ModelIDs in expression:\n")
print(setdiff(tmp$ModelID, expr_sub$ModelID))

if (length(setdiff(tmp$ModelID, expr_sub$ModelID)) > 0) {
    stop("Some ModelIDs were not found in expression file")
}

# =========================================================
# 5. KEEP DEFAULT RNA-seq ENTRY
# =========================================================
if ("IsDefaultEntryForModel" %in% colnames(expr_sub)) {
    expr_sub <- expr_sub[IsDefaultEntryForModel == "Yes"]
}

expr_sub <- unique(expr_sub, by = "ModelID")

# add correct target names and order
expr_sub <- merge(
    expr_sub,
    tmp[, .(
        ModelID,
        CellLineName_final = CellLineName_target,
        StrippedCellLineName,
        order_id
    )],
    by = "ModelID",
    all.x = TRUE
)

setorder(expr_sub, order_id)

cat("\nFinal selected expression rows:\n")
print(expr_sub[, .(ModelID, CellLineName_final, StrippedCellLineName, order_id)])

# =========================================================
# 6. REMOVE METADATA COLUMNS
# =========================================================
meta_cols <- intersect(
    c(
        "V1",
        "SequencingID",
        "ModelID",
        "IsDefaultEntryForModel",
        "ModelConditionID",
        "IsDefaultEntryForMC",
        "CellLineName_final",
        "StrippedCellLineName",
        "order_id"
    ),
    colnames(expr_sub)
)

expr_only <- expr_sub[, !meta_cols, with = FALSE]

# =========================================================
# 7. TRANSPOSE
# Result:
# rows = genes
# columns = selected cell lines
# =========================================================
expr_t <- transpose(expr_only, keep.names = "gene")

colnames(expr_t) <- c("gene", expr_sub$CellLineName_final)

# clean gene names:
# "MICA (100507436)" -> "MICA"
expr_t[, gene := sub(" \\(.*\\)$", "", gene)]

# =========================================================
# 8. FORCE COLUMN ORDER
# =========================================================
final_order <- target_order$CellLineName
existing_order <- final_order[final_order %in% colnames(expr_t)]

expr_t <- expr_t[, c("gene", existing_order), with = FALSE]

# =========================================================
# 9. SAVE FULL RNA-seq MATRIX
# =========================================================
fwrite(expr_t, "DepMap_selected_cell_lines_RNAseq.csv")

# =========================================================
# 10. SAVE MAPPING TABLE
# =========================================================
mapping_out <- tmp[, .(
    CellLineName = CellLineName_target,
    StrippedCellLineName,
    ModelID,
    DepMap_CellLineName = CellLineName_DepMap,
    CCLEName
)]

fwrite(mapping_out, "DepMap_selected_cell_lines_mapping.csv")

# =========================================================
# 11. OPTIONAL: SAVE ONLY NKG2D LIGANDS
# =========================================================
nkg2d_ligands <- c(
    "MICA",
    "MICB",
    "ULBP1",
    "ULBP2",
    "ULBP3",
    "RAET1E",
    "RAET1G",
    "RAET1L"
)

ligand_out <- expr_t[gene %in% nkg2d_ligands]

fwrite(ligand_out, "DepMap_selected_cell_lines_NKG2DL.csv")

# =========================================================
# 12. CHECK OUTPUT
# =========================================================
cat("\nFinal RNA-seq matrix dimensions:\n")
print(dim(expr_t))

cat("\nFirst rows:\n")
print(head(expr_t))

cat("\nNKG2D ligand matrix:\n")
print(ligand_out)

cat("\nSaved files:\n")
cat("1) DepMap_selected_cell_lines_RNAseq.csv\n")
cat("2) DepMap_selected_cell_lines_mapping.csv\n")
cat("3) DepMap_selected_cell_lines_NKG2DL.csv\n")
