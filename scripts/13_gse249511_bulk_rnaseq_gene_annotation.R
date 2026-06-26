# -----------------------------------------------------------------------------
# GSE249511 transcript-to-gene annotation and count preparation
#
# Original archived script: Obajdin et al. [72]/Script annotation.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(rtracklayer)
    library(dplyr)
    library(tibble)
    library(readr)
    library(stringr)
    library(tximport)
})

files <- list.files(
    path = ".",
    pattern = "quant.sf.gz$",
    full.names = TRUE
)

names(files) <- basename(files) %>%
    str_remove("_quant.sf.gz$") %>%
    str_remove("^GSM[0-9]+_") %>%
    str_replace_all("-", "_")

example <- read_tsv(files[1], show_col_types = FALSE)

gtf <- rtracklayer::import("gencode.v49.basic.annotation.gtf")

tx2gene_gtf <- as.data.frame(gtf) %>%
    filter(type == "transcript") %>%
    dplyr::select(transcript_id, gene_name) %>%
    filter(!is.na(transcript_id), !is.na(gene_name)) %>%
    distinct() %>%
    mutate(
        transcript_id_clean = sub("\\..*$", "", transcript_id)
    )

tx2gene_final <- tibble(
    TXNAME = example$Name,
    transcript_id_clean = sub("\\..*$", "", example$Name)
) %>%
    left_join(
        tx2gene_gtf,
        by = "transcript_id_clean"
    ) %>%
    filter(!is.na(gene_name)) %>%
    dplyr::select(
        TXNAME,
        GENEID = gene_name
    )

cat("Total transcripts in quant.sf:\n")
print(length(example$Name))

cat("Mapped transcripts using GTF:\n")
print(nrow(tx2gene_final))

txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene_final,
    countsFromAbundance = "no"
)

counts_gene <- as.data.frame(txi$counts) %>%
    rownames_to_column("SYMBOL")

write_csv(counts_gene, "GSE249511_gene_NumReads_from_raw_quant_GENCODE.csv")
