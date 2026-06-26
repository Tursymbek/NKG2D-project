# -----------------------------------------------------------------------------
# ENCODE TF ChIP-seq promoter overlap heatmaps
#
# Original archived script: after Taiji TF ChIP Seq/NKG2DL_overlap_results/Script with clusters.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

# =========================================================
# TF ChIP-seq overlap with NKG2D ligand promoters
# MCL-grouped heatmap using string_MCL_clusters.tsv
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

chip_dir <- "."
gtf_file <- "gencode.v49.basic.annotation.gtf"
meta_file <- "metadata.tsv"
tf_list_file <- "TF_list_High_enriched_only.txt"

# Use this STRING file
cluster_file <- "string_MCL_clusters.tsv"

out_dir <- file.path(chip_dir, "NKG2DL_overlap_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 1B. Load STRING/MCL clusters and use manual cluster colors
# Legend shows only cluster numbers
# =========================================================

string_clusters_raw <- fread(cluster_file)

string_clusters <- string_clusters_raw %>%
    transmute(
        TF = toupper(trimws(`protein name`)),
        cluster_number = as.numeric(`cluster number`),
        Cluster = as.character(`cluster number`)
    ) %>%
    filter(!is.na(TF), TF != "") %>%
    distinct(TF, .keep_all = TRUE)

# Manual Fill colors from your table
manual_cluster_colors <- c(
    "1"  = "#FF999F",
    "2"  = "#FFD3C1",
    "3"  = "#E0C199",
    "4"  = "#FFF5C1",
    "5"  = "#D7E099",
    "6"  = "#E5FFC1",
    "7"  = "#AFE099",
    "8"  = "#C3FFC1",
    "9"  = "#99E0AB",
    "10" = "#C1FFE3",
    "11" = "#C9C8FF",
    "12" = "#D7C8FF",
    "13" = "#E5C8FF",
    "14" = "#F2C8FF",
    "15" = "#F9C8FF",
    "16" = "#F9C8E3",
    "17" = "#F9C8D5"
)

cluster_levels <- string_clusters %>%
    distinct(Cluster, cluster_number) %>%
    arrange(cluster_number) %>%
    pull(Cluster)

cluster_palette <- manual_cluster_colors[cluster_levels]

cluster_palette <- c(
    cluster_palette
)


# =========================================================
# 2. Target TFs and NKG2D ligands
# =========================================================

high_tfs <- fread(tf_list_file, header = FALSE)$V1
high_tfs <- toupper(trimws(high_tfs))

encode_available_tfs <- c(
    "CEBPB", "ESR1", "ATF4", "ATF2",
    "FOS", "JUN", "RUNX1", "CREB1",
    "SP1", "KLF4", "DDIT3", "JUNB",
    "PBX1", "NFE2L2", "CREB3",
    "CREB3L1", "ATF3", "JUND",
    "FOSL2", "ATF1", "CREM", "MAFG",
    "TBP", "GATA3", "FOSB", "FOSL1",
    "ARNTL", "ATF7", "CEBPD",
    "TFAP2A", "GTF2B", "ETS1",
    "JDP2", "MAFK", "NFE2",
    "TFAP2C", "ATF6", "MAFF",
    "CXXC5", "ARNT2", "MECOM",
    "ELK1", "KLF6", "NFE2L1",
    "BACH1", "HOXA9", "STAT6",
    "CLOCK", "NPAS2", "MYBL2",
    "ETS2", "PBX3"
)

target_tfs <- intersect(high_tfs, encode_available_tfs)

cat("High-enriched TFs total:", length(high_tfs), "\n")
cat("ENCODE available TFs:", length(encode_available_tfs), "\n")
cat("High-enriched TFs available in ENCODE:", length(target_tfs), "\n")
print(target_tfs)

fwrite(
    data.table(TF = target_tfs),
    file.path(out_dir, "TF_list_High_enriched_ENCODE_available.txt"),
    col.names = FALSE
)

nkg2dl_genes <- c(
    "MICA", "MICB",
    "ULBP1", "ULBP2", "ULBP3",
    "RAET1E", "RAET1G", "RAET1L"
)

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
        file_accession = as.character(`File accession`),
        
        TF_clean = toupper(as.character(`Experiment target`)),
        TF_clean = str_replace(TF_clean, "-HUMAN$", ""),
        TF_clean = str_replace(TF_clean, "\\s+.*$", ""),
        
        cell_line = as.character(`Biosample term name`),
        output_type = as.character(`Output type`),
        assay = as.character(`Assay`),
        file_format = as.character(`File format`)
    ) %>%
    select(
        file_accession,
        TF_clean,
        cell_line,
        assay,
        output_type,
        file_format
    )

meta2 <- meta2 %>%
    filter(
        TF_clean %in% target_tfs,
        assay == "TF ChIP-seq" | assay == "ChIP-seq",
        file_format %in% c("bed", "bed narrowPeak", "bed broadPeak", "bed gappedPeak") |
            grepl("peak", tolower(output_type))
    )

cat("Metadata files after TF/assay/peak filtering:", nrow(meta2), "\n")
print(table(meta2$TF_clean))

# =========================================================
# 5. Find local BED/GZ peak files
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

cat("Found local BED files:", length(bed_files), "\n")

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
    stop("No files matched target TFs. Check metadata.tsv and downloaded BED files.")
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
# 10. Matrix TF x NKG2DL: total peaks
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
# 11. Normalized matrix: peaks per file
# =========================================================

file_counts <- file_anno_use %>%
    mutate(TF = toupper(as.character(TF_clean))) %>%
    count(TF, name = "n_files_total")

summary_norm <- summary_tf_gene %>%
    mutate(TF = toupper(as.character(TF))) %>%
    left_join(file_counts, by = "TF") %>%
    mutate(
        peaks_per_file = total_overlapping_peaks / n_files_total
    )

write.csv(
    summary_norm,
    file.path(out_dir, "SUMMARY_TF_by_NKG2DL_gene_normalized.csv"),
    row.names = FALSE
)

norm_matrix <- summary_norm %>%
    select(TF, gene_name, peaks_per_file) %>%
    pivot_wider(
        names_from = gene_name,
        values_from = peaks_per_file,
        values_fill = 0
    )

missing_tfs <- setdiff(target_tfs, norm_matrix$TF)

if (length(missing_tfs) > 0) {
    tmp <- data.frame(TF = missing_tfs)
    for (g in nkg2dl_genes) tmp[[g]] <- 0
    norm_matrix <- bind_rows(norm_matrix, tmp)
}

for (g in nkg2dl_genes) {
    if (!g %in% colnames(norm_matrix)) {
        norm_matrix[[g]] <- 0
    }
}

norm_matrix <- norm_matrix %>%
    mutate(TF = factor(TF, levels = target_tfs)) %>%
    arrange(TF)

write.csv(
    norm_matrix,
    file.path(out_dir, "MATRIX_TF_x_NKG2DL_peaks_per_file.csv"),
    row.names = FALSE
)

# =========================================================
# 12. Helper functions for STRING-cluster-grouped heatmaps
# =========================================================

make_mcl_order <- function(mat_input) {
    
    row_info <- data.frame(
        TF = rownames(mat_input),
        stringsAsFactors = FALSE
    ) %>%
        left_join(
            string_clusters %>%
                select(TF, Cluster, cluster_number),
            by = "TF"
        ) %>%
        mutate(
            Cluster = ifelse(is.na(Cluster), "Not in STRING", Cluster),
            cluster_number = ifelse(is.na(cluster_number), 999, cluster_number)
        ) %>%
        arrange(cluster_number, TF)
    
    ordered_tfs <- c()
    
    cluster_order <- unique(row_info$Cluster)
    
    for (cl in cluster_order) {
        
        tfs_in_cluster <- row_info$TF[row_info$Cluster == cl]
        
        if (length(tfs_in_cluster) > 2) {
            
            submat <- mat_input[tfs_in_cluster, , drop = FALSE]
            row_sd <- apply(submat, 1, sd, na.rm = TRUE)
            
            if (sum(row_sd > 0) > 1) {
                
                hc <- hclust(
                    dist(submat, method = "euclidean"),
                    method = "complete"
                )
                
                ordered_tfs <- c(ordered_tfs, tfs_in_cluster[hc$order])
                
            } else {
                ordered_tfs <- c(ordered_tfs, tfs_in_cluster)
            }
            
        } else {
            ordered_tfs <- c(ordered_tfs, tfs_in_cluster)
        }
    }
    
    row_info_ordered <- row_info[
        match(ordered_tfs, row_info$TF),
        ,
        drop = FALSE
    ]
    
    row_annot <- data.frame(
        Cluster = row_info_ordered$Cluster,
        row.names = row_info_ordered$TF,
        stringsAsFactors = FALSE
    )
    
    row_annot$Cluster <- factor(
        row_annot$Cluster,
        levels = names(cluster_palette)
    )
    
    cluster_runs <- rle(as.character(row_annot$Cluster))
    gaps_row <- cumsum(cluster_runs$lengths)
    
    if (length(gaps_row) > 1) {
        gaps_row <- gaps_row[-length(gaps_row)]
    } else {
        gaps_row <- NULL
    }
    
    return(list(
        ordered_tfs = ordered_tfs,
        row_annot = row_annot,
        gaps_row = gaps_row
    ))
}


plot_mcl_heatmap <- function(
        mat_input,
        output_file,
        main_title,
        color_palette,
        number_format = "%.1f",
        legend_breaks = NA,
        legend_labels = NA
) {
    
    order_obj <- make_mcl_order(mat_input)
    
    mat_ordered <- mat_input[
        order_obj$ordered_tfs,
        ,
        drop = FALSE
    ]
    
    annotation_colors <- list(
        Cluster = cluster_palette
    )
    
    png(
        file.path(out_dir, output_file),
        width = 3000,
        height = 3600,
        res = 300
    )
    
    pheatmap(
        mat_ordered,
        
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        
        gaps_row = order_obj$gaps_row,
        
        annotation_row = order_obj$row_annot,
        annotation_colors = annotation_colors,
        
        color = color_palette,
        border_color = "grey85",
        
        fontsize = 15,
        fontsize_row = 13,
        fontsize_col = 18,
        angle_col = 45,
        
        display_numbers = TRUE,
        number_format = number_format,
        fontsize_number = 12,
        number_color = "black",
        
        legend_breaks = legend_breaks,
        legend_labels = legend_labels,
        
        main = main_title
    )
    
    dev.off()
}


# =========================================================
# 13. Prepare matrices for plotting
# =========================================================

mat_total <- as.data.frame(binding_matrix)
rownames(mat_total) <- as.character(mat_total$TF)
mat_total$TF <- NULL
mat_total <- mat_total[, nkg2dl_genes, drop = FALSE]
mat_total <- as.matrix(mat_total)
colnames(mat_total) <- nkg2dl_labels

mat_log <- log2(mat_total + 1)
colnames(mat_log) <- nkg2dl_labels

mat_bin <- ifelse(mat_total > 0, 1, 0)
colnames(mat_bin) <- nkg2dl_labels

mat_norm <- as.data.frame(norm_matrix)
rownames(mat_norm) <- as.character(mat_norm$TF)
mat_norm$TF <- NULL
mat_norm <- mat_norm[, nkg2dl_genes, drop = FALSE]
mat_norm <- as.matrix(mat_norm)
colnames(mat_norm) <- nkg2dl_labels


# =========================================================
# 14. Plot STRING-cluster-grouped heatmaps
# =========================================================

plot_mcl_heatmap(
    mat_input = mat_log,
    output_file = "FINAL_heatmap_TF_NKG2DL_log2_STRING_CLUSTER_GROUPED.png",
    main_title = "High-enriched TF ChIP-seq binding to NKG2D ligand promoters\nlog2(overlapping peaks + 1)",
    color_palette = colorRampPalette(
        c("white", "#FEE08B", "#F46D43", "#A50026")
    )(100),
    number_format = "%.1f"
)

plot_mcl_heatmap(
    mat_input = mat_bin,
    output_file = "FINAL_heatmap_TF_NKG2DL_binary_STRING_CLUSTER_GROUPED.png",
    main_title = "High-enriched TF ChIP-seq binding to NKG2D ligand promoters",
    color_palette = c("white", "#E63946"),
    number_format = "%.0f",
    legend_breaks = c(0, 1),
    legend_labels = c("No binding", "Binding")
)

plot_mcl_heatmap(
    mat_input = mat_norm,
    output_file = "FINAL_heatmap_TF_NKG2DL_normalized_peaks_per_file_STRING_CLUSTER_GROUPED.png",
    main_title = "High-enriched TF ChIP-seq binding to NKG2D ligand promoters\nnormalized peaks per file",
    color_palette = colorRampPalette(
        c("white", "#FEE08B", "#F46D43", "#A50026")
    )(100),
    number_format = "%.1f"
)

cat("\n✅ DONE\n")
cat("Saved outputs to:", out_dir, "\n")
