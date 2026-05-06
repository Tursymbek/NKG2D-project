# =========================================================
# TF ChIP-seq TOP30 overlap with NKG2D ligand promoters
# Works with: D:/NKG2D project/TF ChIP-Seq TOP30
# =========================================================

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(GenomicRanges)
    library(rtracklayer)
    library(stringr)
    library(pheatmap)
})

# =========================================================
# 1. Paths
# =========================================================

chip_dir <- "D:/NKG2D project/TF ChIP-Seq TOP30"

gtf_file <- "D:/NKG2D project/TF ChIP-Seq TOP30/gencode.v49.basic.annotation.gtf"
meta_file <- file.path(chip_dir, "metadata.tsv")

out_dir <- file.path(chip_dir, "NKG2DL_overlap_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. Target TFs and NKG2D ligands
# =========================================================

target_tfs <- c(
    "SP1", "TP53", "CEBPB", "CREB1", "SRF", "RARA",
    "RXRA", "STAT3", "RUNX1", "SPI1", "TAL1", "TCF12",
    "CEBPA", "TCF4", "POU5F1", "TCF3", "MYC", "RXRB",
    "GATA4", "E2F1", "ATF4", "FOXO1", "SREBF1"
)

nkg2dl_genes <- c(
    "MICA", "MICB",
    "ULBP1", "ULBP2", "ULBP3",
    "RAET1E", "RAET1G", "RAET1L"
)

# labels shown on heatmaps
nkg2dl_labels <- c(
    "MICA", "MICB",
    "ULBP1", "ULBP2", "ULBP3",
    "ULBP4", "ULBP5", "ULBP6"
)

# =========================================================
# 3. Load GTF and create promoters
# =========================================================

cat("Loading GTF...\n")

gtf <- import(gtf_file)

genes <- gtf[gtf$type == "gene"]
genes <- genes[genes$gene_name %in% nkg2dl_genes]

if (length(genes) == 0) {
    stop("No NKG2D ligand genes found in GTF. Check gene names or GTF file.")
}

promoters_gr <- promoters(
    genes,
    upstream = 2000,
    downstream = 2000
)

promoters_gr$gene_name <- genes$gene_name
promoters_gr$gene_id <- genes$gene_id
promoters_gr$score <- 0

promoter_bed <- file.path(out_dir, "NKG2DL_promoters_TSS_plusminus_2kb.bed")
export(promoters_gr, promoter_bed, format = "BED")

cat("Promoters saved to:", promoter_bed, "\n")

# =========================================================
# 4. Load metadata
# =========================================================

meta <- fread(meta_file)

cat("Metadata columns:\n")
print(colnames(meta))

meta2 <- meta %>%
    mutate(
        file_accession = as.character(File),
        TF_clean = toupper(as.character(TF)),
        TF_clean = str_replace(TF_clean, "-HUMAN$", ""),
        TF_clean = str_replace(TF_clean, "\\s+.*$", ""),
        cell_line = as.character(CellLine),
        output_type = as.character(OutputType),
        assay = as.character(Assay)
    ) %>%
    select(file_accession, TF_clean, cell_line, assay, output_type)

# =========================================================
# 5. Find local BED.GZ files
# =========================================================

bed_files <- list.files(
    chip_dir,
    pattern = "\\.bed\\.gz$|\\.bed$|\\.narrowPeak\\.gz$|\\.broadPeak\\.gz$",
    full.names = TRUE,
    recursive = FALSE
)

if (length(bed_files) == 0) {
    stop("No BED / BED.GZ peak files found in chip_dir.")
}

cat("Found BED files:", length(bed_files), "\n")

get_accession <- function(x) {
    str_extract(basename(x), "ENCFF[0-9A-Z]+")
}

file_anno <- data.frame(
    file = bed_files,
    file_name = basename(bed_files),
    file_accession = get_accession(bed_files),
    stringsAsFactors = FALSE
) %>%
    left_join(meta2, by = "file_accession")

write.csv(
    file_anno,
    file.path(out_dir, "local_ChIPseq_file_annotation.csv"),
    row.names = FALSE
)

cat("Files with TF annotation:", sum(!is.na(file_anno$TF_clean)), "\n")

# =========================================================
# 6. Keep only selected TFs
# =========================================================

file_anno_use <- file_anno %>%
    filter(
        !is.na(TF_clean),
        TF_clean %in% target_tfs
    )

cat("Files matching target TFs:", nrow(file_anno_use), "\n")
print(table(file_anno_use$TF_clean))

if (nrow(file_anno_use) == 0) {
    stop("No files matched target TFs. Check metadata.tsv.")
}

# =========================================================
# 7. Function to read BED peak files
# =========================================================

read_peak_file <- function(f) {
    dt <- fread(f, header = FALSE)
    
    if (ncol(dt) < 3) {
        stop(paste("File has <3 columns:", f))
    }
    
    dt <- dt[, 1:3]
    colnames(dt) <- c("chr", "start", "end")
    
    dt <- dt %>%
        filter(
            !is.na(chr),
            !is.na(start),
            !is.na(end),
            str_detect(chr, "^chr"),
            end > start
        )
    
    GRanges(
        seqnames = dt$chr,
        ranges = IRanges(start = dt$start + 1, end = dt$end)
    )
}

# =========================================================
# 8. Overlap analysis
# =========================================================

cat("Running overlap analysis...\n")

all_hits <- list()

for (i in seq_len(nrow(file_anno_use))) {
    
    f <- file_anno_use$file[i]
    
    cat("Processing:", i, "/", nrow(file_anno_use), basename(f), "\n")
    
    peaks_gr <- tryCatch(
        read_peak_file(f),
        error = function(e) {
            message("Failed to read: ", f)
            message(e$message)
            return(NULL)
        }
    )
    
    if (is.null(peaks_gr)) next
    
    hits <- findOverlaps(promoters_gr, peaks_gr)
    
    if (length(hits) == 0) {
        
        tmp <- data.frame(
            file_name = basename(f),
            file_accession = file_anno_use$file_accession[i],
            TF = file_anno_use$TF_clean[i],
            cell_line = file_anno_use$cell_line[i],
            gene_name = NA_character_,
            n_overlapping_peaks = 0
        )
        
    } else {
        
        hit_df <- data.frame(
            gene_name = promoters_gr$gene_name[queryHits(hits)],
            peak_id = subjectHits(hits)
        )
        
        tmp <- hit_df %>%
            group_by(gene_name) %>%
            summarise(
                n_overlapping_peaks = n_distinct(peak_id),
                .groups = "drop"
            ) %>%
            mutate(
                file_name = basename(f),
                file_accession = file_anno_use$file_accession[i],
                TF = file_anno_use$TF_clean[i],
                cell_line = file_anno_use$cell_line[i]
            ) %>%
            select(
                file_name,
                file_accession,
                TF,
                cell_line,
                gene_name,
                n_overlapping_peaks
            )
    }
    
    all_hits[[length(all_hits) + 1]] <- tmp
}

overlap_results <- bind_rows(all_hits)

write.csv(
    overlap_results,
    file.path(out_dir, "TF_ChIPseq_NKG2DL_promoter_overlap_long.csv"),
    row.names = FALSE
)

# =========================================================
# 9. Summary table
# =========================================================

summary_tf_gene <- overlap_results %>%
    filter(!is.na(gene_name)) %>%
    group_by(TF, gene_name) %>%
    summarise(
        n_files_with_binding = sum(n_overlapping_peaks > 0),
        total_overlapping_peaks = sum(n_overlapping_peaks, na.rm = TRUE),
        cell_lines = paste(sort(unique(na.omit(cell_line))), collapse = "; "),
        file_accessions = paste(sort(unique(na.omit(file_accession))), collapse = "; "),
        .groups = "drop"
    ) %>%
    arrange(TF, gene_name)

write.csv(
    summary_tf_gene,
    file.path(out_dir, "SUMMARY_TF_by_NKG2DL_gene.csv"),
    row.names = FALSE
)

# =========================================================
# 10. Matrix TF x NKG2DL
# =========================================================

binding_matrix <- summary_tf_gene %>%
    select(TF, gene_name, total_overlapping_peaks) %>%
    pivot_wider(
        names_from = gene_name,
        values_from = total_overlapping_peaks,
        values_fill = 0
    )

missing_tfs <- setdiff(target_tfs, binding_matrix$TF)

if (length(missing_tfs) > 0) {
    tmp <- data.frame(TF = missing_tfs)
    for (g in nkg2dl_genes) tmp[[g]] <- 0
    binding_matrix <- bind_rows(binding_matrix, tmp)
}

for (g in nkg2dl_genes) {
    if (!g %in% colnames(binding_matrix)) {
        binding_matrix[[g]] <- 0
    }
}

binding_matrix <- binding_matrix %>%
    mutate(TF = factor(TF, levels = target_tfs)) %>%
    arrange(TF)

write.csv(
    binding_matrix,
    file.path(out_dir, "MATRIX_TF_x_NKG2DL_total_peaks.csv"),
    row.names = FALSE
)

# =========================================================
# 11. Heatmap: total peaks log2
# =========================================================

mat <- as.data.frame(binding_matrix)
rownames(mat) <- as.character(mat$TF)
mat$TF <- NULL
mat <- mat[, nkg2dl_genes, drop = FALSE]

mat_log <- log2(as.matrix(mat) + 1)
colnames(mat_log) <- nkg2dl_labels

png(
    file.path(out_dir, "FINAL_heatmap_TF_NKG2DL_TOP30_log2.png"),
    width = 3000,
    height = 2200,
    res = 300
)

pheatmap(
    mat_log,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "#FEE08B", "#F46D43", "#A50026"))(100),
    border_color = "grey80",
    fontsize = 18,
    fontsize_row = 18,
    fontsize_col = 18,
    angle_col = 45,
    main = "TOP TF ChIP-seq binding to NKG2D ligand promoters\nlog2(overlapping peaks + 1)"
)

dev.off()

# =========================================================
# 12. Binary heatmap
# =========================================================

mat_bin <- ifelse(mat > 0, 1, 0)
colnames(mat_bin) <- nkg2dl_labels

png(
    file.path(out_dir, "FINAL_heatmap_TF_NKG2DL_TOP30_binary.png"),
    width = 3000,
    height = 2200,
    res = 300
)

pheatmap(
    mat_bin,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = c("white", "#E63946"),
    border_color = "grey80",
    fontsize = 18,
    fontsize_row = 18,
    fontsize_col = 18,
    angle_col = 45,
    legend_breaks = c(0, 1),
    legend_labels = c("No binding", "Binding"),
    main = "TOP TF ChIP-seq binding to NKG2D ligand promoters"
)

dev.off()

# =========================================================
# 13. Normalized matrix: peaks per file
# =========================================================

# Используем ту же колонку TF, что уже есть в summary_tf_gene
file_counts <- file_anno_use %>%
    dplyr::mutate(TF = toupper(as.character(TF_clean))) %>%
    dplyr::count(TF, name = "n_files_total")

# Проверка
print(file_counts)

summary_norm <- summary_tf_gene %>%
    dplyr::mutate(TF = toupper(as.character(TF))) %>%
    dplyr::left_join(file_counts, by = "TF") %>%
    dplyr::mutate(
        peaks_per_file = total_overlapping_peaks / n_files_total
    )

write.csv(
    summary_norm,
    file.path(out_dir, "SUMMARY_TF_by_NKG2DL_gene_normalized.csv"),
    row.names = FALSE
)

# =========================================================
# Matrix
# =========================================================

norm_matrix <- summary_norm %>%
    dplyr::select(TF, gene_name, peaks_per_file) %>%
    tidyr::pivot_wider(
        names_from = gene_name,
        values_from = peaks_per_file,
        values_fill = 0
    )

# добавить отсутствующие TF
missing_tfs <- setdiff(target_tfs, norm_matrix$TF)

if (length(missing_tfs) > 0) {
    tmp <- data.frame(TF = missing_tfs)
    for (g in nkg2dl_genes) tmp[[g]] <- 0
    norm_matrix <- dplyr::bind_rows(norm_matrix, tmp)
}

# добавить отсутствующие гены
for (g in nkg2dl_genes) {
    if (!g %in% colnames(norm_matrix)) {
        norm_matrix[[g]] <- 0
    }
}

norm_matrix <- norm_matrix %>%
    dplyr::mutate(TF = factor(TF, levels = target_tfs)) %>%
    dplyr::arrange(TF)

write.csv(
    norm_matrix,
    file.path(out_dir, "MATRIX_TF_x_NKG2DL_peaks_per_file.csv"),
    row.names = FALSE
)

# =========================================================
# 14. Heatmap (normalized)
# =========================================================

mat_norm <- as.data.frame(norm_matrix)
rownames(mat_norm) <- as.character(mat_norm$TF)
mat_norm$TF <- NULL
mat_norm <- mat_norm[, nkg2dl_genes, drop = FALSE]

# правильные подписи
colnames(mat_norm) <- c(
    "MICA", "MICB",
    "ULBP1", "ULBP2", "ULBP3",
    "ULBP4", "ULBP5", "ULBP6"
)

png(
    file.path(out_dir, "FINAL_heatmap_TF_NKG2DL_TOP30_normalized_peaks_per_file.png"),
    width = 3000,
    height = 2200,
    res = 300
)

pheatmap(
    as.matrix(mat_norm),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "#FEE08B", "#F46D43", "#A50026"))(100),
    border_color = "grey80",
    fontsize = 18,
    fontsize_row = 18,
    fontsize_col = 18,
    angle_col = 45,
    main = "TOP TF ChIP-seq binding to NKG2D ligand promoters\nnormalized peaks per file"
)

dev.off()

cat("\n✅ NORMALIZED DONE\n")