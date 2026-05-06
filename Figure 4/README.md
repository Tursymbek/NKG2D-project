# TF ChIP-seq Promoter Binding Analysis of NKG2D Ligands

## Overview

This repository contains an R-based pipeline for the analysis of transcription factor (TF) binding at promoter regions of NKG2D ligand genes using ENCODE ChIP-seq data. The analysis quantifies TF binding by computing overlaps between TF peak regions and promoter intervals defined as transcription start sites (TSS) ± 2 kb.

The primary objective of this workflow is to identify transcription factors that potentially regulate the expression of NKG2D ligands, including MICA, MICB, and members of the ULBP family (ULBP1–6, corresponding to RAET1E, RAET1G, and RAET1L in GENCODE annotation).

---

## Selection of Transcription Factors

The initial list of candidate transcription factors was derived from Taiji analysis and included the following TOP30 TFs:

TP53, STAT3, SP1, RUNX1, RXRA, CEBPB, RARA, CREB1, SPI1, SRF, TAL1, TCF12, CEBPA, TCF4, POU5F1, RUNX3, FLI1, TCF3, MYC, KLF5, RUNX2, RXRG, RXRB, GATA4, FOXO3, ASCL1, E2F1, ATF4, FOXO1, SREBF1.

However, not all transcription factors had available or suitable ChIP-seq datasets in the ENCODE database under the applied filtering criteria (human, cancer cell lines, processed peak files).

Therefore, the final analysis was restricted to transcription factors for which ENCODE ChIP-seq peak files were available and successfully downloaded:

SP1, TP53, CEBPB, CREB1, SRF, RARA, RXRA, STAT3, RUNX1, SPI1, TAL1, TCF12, CEBPA, TCF4, POU5F1, TCF3, MYC, RXRB, GATA4, E2F1, ATF4, FOXO1, SREBF1.

Transcription factors excluded from the analysis due to lack of suitable ENCODE datasets included:

RUNX3, FLI1, KLF5, RUNX2, RXRG, FOXO3, ASCL1.

This filtering step ensures that all analyzed TFs are supported by experimentally validated ChIP-seq data within a consistent dataset framework.
---

## Data Source and Preprocessing

The analysis is based on processed ChIP-seq peak files downloaded from the ENCODE portal. The dataset was restricted to human TF ChIP-seq experiments performed in cancer cell lines, and only processed peak files in `.bed.gz` format were retained for downstream analysis. Other file types such as BAM, bigWig, and bigBed were excluded, as they are not required for genomic overlap analysis.

The list of downloaded files was stored in `LIST_of_files.txt` and manually filtered to include only the files that were actually used.

---

## Metadata Processing

Metadata for each ChIP-seq file were retrieved through the ENCODE REST API using file accession identifiers. JSON responses were parsed to extract relevant fields, including:

- Transcription factor name  
- Cell line  
- Assay type  
- Laboratory source  
- Output type  

These data were consolidated into a single file:

```

metadata.tsv

```

Matching between metadata and peak files was performed using ENCODE file accession IDs extracted from file names.

---

## Cell Lines

The final dataset includes experiments from the following cancer cell lines:

- K562  
- HepG2  
- A549  
- MCF-7  
- HeLa-S3  
- Ishikawa  
- HCT116  
- SK-N-SH  
- HL-60  
- NB4  

Due to uneven representation across transcription factors, normalization by the number of files per TF was implemented.

---

## Gene Annotation

Gene annotation was performed using:

```

GENCODE v49 (GRCh38)

```

Promoter regions were defined as:

```

TSS ± 2000 bp

````

The analyzed genes include:

- MICA, MICB  
- ULBP1, ULBP2, ULBP3  
- RAET1E (ULBP4)  
- RAET1G (ULBP5)  
- RAET1L (ULBP6)  

Minor differences between recent GENCODE versions are not expected to significantly affect results.

---

## Analysis Workflow

1. **Load annotation**  
   GTF file is imported and gene coordinates are extracted.

2. **Generate promoter regions**  
   Promoters are defined using ±2 kb around TSS.

3. **Load ChIP-seq peaks**  
   Peak files (`.bed.gz`, `.narrowPeak`, `.broadPeak`) are read and converted to genomic ranges.

4. **Overlap analysis**  
   TF peaks are intersected with promoter regions using:

   ```r
   findOverlaps(promoters_gr, peaks_gr)
   ```

5. **Aggregation**
   Overlapping peaks are counted for each TF–gene pair across all files and cell lines.

---

## Output Files 

### Long-format results

```
TF_ChIPseq_NKG2DL_promoter_overlap_long.csv
```

Contains:

* TF
* Gene
* Cell line
* Number of overlapping peaks

---

### Summary table

```
SUMMARY_TF_by_NKG2DL_gene.csv
```

Includes:

* Total overlapping peaks
* Number of files with binding
* Cell lines
* File accessions

---

### TF × Gene matrix

```
MATRIX_TF_x_NKG2DL_total_peaks.csv
```

---

### Normalized matrix

```
MATRIX_TF_x_NKG2DL_peaks_per_file.csv
```

Normalization is performed as:

```r
peaks_per_file = total_overlapping_peaks / n_files_total
```

This step corrects for unequal numbers of ChIP-seq experiments per transcription factor.

---

## Visualization

Heatmaps are generated using the `pheatmap` package and include:

* Log2-transformed peak counts
* Binary binding (presence/absence)
* Normalized peaks per file

Gene labels are standardized to ULBP1–6 for clarity.

---

## Environment

```
R version 4.5.3 (2026-03-11)
Platform: Windows 11 x64
```

### Key packages

* data.table
* dplyr
* tidyr
* stringr
* GenomicRanges
* IRanges
* rtracklayer
* pheatmap

---

## Reproducibility

To reproduce the analysis:

1. Download ENCODE `.bed.gz` peak files
2. Generate `metadata.tsv` via ENCODE API
3. Download GENCODE v49 annotation
4. Place all files in the working directory
5. Run:

```r
source("TF_ChIPseq_NKG2DL_promoter_overlap.R")
```

The script will automatically detect input files, perform overlap analysis, and generate all output tables and figures.

---

## Notes

* Only processed peak files (`.bed.gz`) were used
* Results depend on dataset availability and TF coverage
* Promoter-based analysis is robust to minor annotation differences

---
