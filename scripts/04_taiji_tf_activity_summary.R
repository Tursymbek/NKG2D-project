# -----------------------------------------------------------------------------
# Taiji transcription-factor PageRank summary
#
# Original archived script: Taiji/Script.txt
# Repository note: converted from .txt to .R and lightly path-normalized.
# Raw public datasets are not bundled unless small derived outputs were present
# in the deposited archive. See README.md and data/external/README.md.
# -----------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================================
# 1. LOAD DATA
# =========================================================
ranks <- fread("GeneRanks.tsv", check.names = FALSE)
setnames(ranks, 1, "TF")

# =========================================================
# 2. TAKE FIRST HIGH AND FIRST LOW COLUMN
# =========================================================
high_idx <- which(colnames(ranks) == "High")[1]
low_idx  <- which(colnames(ranks) == "Low")[1]

ranks_simple <- data.table(
    TF   = ranks[[1]],
    High = as.numeric(ranks[[high_idx]]),
    Low  = as.numeric(ranks[[low_idx]])
)

# =========================================================
# 3. HIGH vs LOW
# =========================================================
eps <- 1e-10

res <- ranks_simple %>%
    mutate(
        log2FC_High_vs_Low = log2((High + eps) / (Low + eps)),
        diff_High_minus_Low = High - Low,
        Regulation = case_when(
            log2FC_High_vs_Low > 0 ~ "Higher in High / Lower in Low",
            log2FC_High_vs_Low < 0 ~ "Higher in Low / Lower in High",
            TRUE ~ "No change"
        )
    ) %>%
    arrange(desc(log2FC_High_vs_Low))

# =========================================================
# 4. SAVE FULL RESULTS
# =========================================================
fwrite(res, "TF_High_vs_Low_simple.tsv", sep = "\t")

high_enriched <- res %>%
    filter(log2FC_High_vs_Low > 0) %>%
    arrange(desc(log2FC_High_vs_Low))

low_enriched <- res %>%
    filter(log2FC_High_vs_Low < 0) %>%
    arrange(log2FC_High_vs_Low)

fwrite(high_enriched, "TFs_increased_in_High_decreased_in_Low.tsv", sep = "\t")
fwrite(low_enriched,  "TFs_increased_in_Low_decreased_in_High.tsv", sep = "\t")

fwrite(head(high_enriched, 30), "TOP30_TFs_High.tsv", sep = "\t")
fwrite(head(low_enriched, 30),  "TOP30_TFs_Low.tsv", sep = "\t")

# =========================================================
# 5. SAVE CLEAN TF LISTS ONLY
# =========================================================
fwrite(
    data.table(TF = high_enriched$TF),
    "TF_list_High_enriched_only.txt",
    col.names = FALSE
)

fwrite(
    data.table(TF = low_enriched$TF),
    "TF_list_Low_enriched_only.txt",
    col.names = FALSE
)

fwrite(
    data.table(TF = head(high_enriched$TF, 30)),
    "TOP30_TF_list_High_only.txt",
    col.names = FALSE
)

fwrite(
    data.table(TF = head(low_enriched$TF, 30)),
    "TOP30_TF_list_Low_only.txt",
    col.names = FALSE
)
