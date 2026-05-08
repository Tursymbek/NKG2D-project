## TCAT / starCAT analysis of speedingCARs single-cell RNA-seq dataset (GSE214231)

This repository contains R scripts used for preprocessing, TCAT/starCAT annotation, NKG2D ligand analysis, UMAP visualization, CAR composition analysis, and cell cycle association analysis of the publicly available single-cell RNA-seq dataset **GSE214231** from the study:

> Castellanos-Rueda R, Di Roberto RB, Biebrich F, Schlatter FS et al.
> *speedingCARs: accelerating the engineering of CAR T cells by signaling domain shuffling and single-cell sequencing.*
> Nature Communications. 2022;13(1):6555. https://doi.org/10.1038/s41467-022-34141-8
 

Dataset source:
GEO accession: **GSE214231**

>Kotliar, D., Curtis, M., Agnew, R., Weinand, K., Nathan, A., Baglaenko, Y., '
>Slowikowski, K., Zhao, Y., Sabeti, P. C., Rao, D. A., & Raychaudhuri, S. (2025). '
>Reproducible single-cell annotation of programs underlying T cell subsets, activation states and functions.
>Nature methods, 22(9), 1964–1980. https://doi.org/10.1038/s41592-025-02793-1


---

# Dataset description

The original study generated a library of approximately 180 unique CAR variants integrated into primary human T cells using CRISPR-Cas9 engineering. Functional pooled screening was performed by co-culture with HER2-positive SKBR3 breast cancer cells followed by:

* single-cell RNA sequencing (scRNA-seq)
* single-cell CAR sequencing (scCAR-seq)

The processed data were downloaded from GEO.

---

# GEO files used

The following files were used from GEO:

## Main GEO accession

* **GSE214231**

## Processed supplementary archive

* `GSE214231_RAW.tar`

## Barcode annotation files

Used for assigning CAR identities to single cells:

* `GSM6601767_D1_2DomLib_barcodes_assigned.csv`
* `GSM6601768_D1_2DomLib_WT_barcodes_assigned.csv`
* `GSM6601769_D2_2DomLib_WT_barcodes_assigned.csv`
* `GSM6601770_D3_2DomLib_1_barcodes_assigned.csv`
* `GSM6601771_D3_2DomLib_2_barcodes_assigned.csv`
* `GSM6601772_D3_2DomLib_WT_barcodes_assigned.csv`

## 10x Genomics matrices

The following sample folders were imported with `Read10X()`:

* `D1_2DomLib`
* `D1_2DomLib_WT`
* `D2_2DomLib_WT`
* `D2_Unstim_28z`
* `D3_2DomLib_1`
* `D3_2DomLib_2`
* `D3_2DomLib_WT`
* `D3_hT`

---

# Analysis overview

The workflow consisted of:

1. Importing and merging 10x Genomics scRNA-seq datasets
2. Assigning CAR variants using barcode annotation tables
3. Quality control filtering
4. T-cell enrichment and tumor-cell exclusion
5. Seurat normalization and dimensional reduction
6. Export to starCAT-compatible format
7. TCAT/starCAT program annotation
8. NKG2D ligand expression analysis
9. UMAP and CAR composition analysis
10. Cell cycle scoring
11. Correlation analysis between TCAT programs and NKG2D ligand expression

---

# Software environment

Analysis was performed in:

* R version 4.5.3
* Windows 11 x64

Main packages:

* Seurat 5.4.0
* SeuratObject 5.3.0
* dplyr
* ggplot2
* uwot
* scatterpie
* Matrix
* tidyr
* readr

---

# NKG2D ligands analyzed

The following NKG2D ligands were included:

| Gene symbol | Alternative name |
| ----------- | ---------------- |
| MICA        | MICA             |
| MICB        | MICB             |
| ULBP1       | ULBP1            |
| ULBP2       | ULBP2            |
| ULBP3       | ULBP3            |
| RAET1E      | ULBP4            |
| RAET1G      | ULBP5            |
| RAET1L      | ULBP6            |

---

# Script descriptions

---

# 1. Data preparation and starCAT analysis

File:
`Data preparation - TCAT.txt` 

## Purpose

This script performs:

* loading all 10x datasets
* CAR barcode assignment
* Seurat preprocessing
* T-cell filtering
* export into starCAT format
* launching starCAT
* importing TCAT program scores back into Seurat

## Main steps

### Dataset import

All GEO sample folders are loaded using:

```r
Read10X()
CreateSeuratObject()
```

### CAR assignment

CAR identities are assigned using barcode annotation CSV files.

### Quality control

Cells were filtered using:

```r
nCount_RNA >= 500
nFeature_RNA >= 200
nFeature_RNA <= 7000
percent.mt <= 15
```

### T-cell enrichment

T cells were enriched using canonical T-cell markers:

* PTPRC
* CD3D
* CD3E
* TRBC1
* TRBC2
* IL7R
* NKG7
* CD4
* CD8A

Tumor-associated cells were removed using:

* EPCAM
* ERBB2
* KRT8
* KRT18
* KRT19

### Seurat preprocessing

The following standard Seurat workflow was applied:

```r
NormalizeData()
FindVariableFeatures()
ScaleData()
RunPCA()
RunUMAP()
FindNeighbors()
FindClusters()
```

### starCAT export

The filtered T-cell matrix was exported as:

* `matrix.mtx.gz`
* `barcodes.tsv.gz`
* `features.tsv.gz`

### starCAT execution

starCAT was launched using:

```bash
starcat --reference "TCAT.V1"
```

### TCAT integration

TCAT scores and dominant programs were re-imported into Seurat metadata.

## Main outputs

* `seurat_T_with_starCAT.rds`
* `metadata_cells_with_starCAT.csv`
* `GSE214231_Tcells.scores.txt`

---

# 2. UMAP and CAR composition analysis

File:
`UMAP, Composition.txt` 

## Purpose

This script generates:

* global UMAP visualization
* CAR-specific UMAPs
* TCAT composition plots
* CAR ranking tables

## Main analyses

### CAR filtering

Only CAR variants with at least:

```r
>= 100 cells
```

were retained.

### UMAP generation

UMAP embeddings were generated using:

* ASA scores
* proliferation scores
* TCAT multinomial labels

via the `uwot` package.

### Composition analysis

TCAT label composition was calculated for each CAR variant.

### CAR ranking

CARs were ranked according to:

* mean ASA score
* mean proliferation score
* dominant TCAT label

## Main outputs

* `UMAP_global_by_CAR.png`
* `UMAP_global_by_TCAT_label.png`
* `TCAT_cluster_composition_by_CAR.png`
* `CAR_ranking_TCAT.csv`
* `UMAP_each_CAR/`

---

# 3. NKG2D ligand score by CAR

File:
`NKG2D ligand SUM score by CAR.txt` 

## Purpose

This script calculates total NKG2D ligand expression for each CAR variant.

## Method

The total ligand score was calculated as:

NKG2DL_{sum}=\sum_{i=1}^{n}Expression_{ligand_i}

using:

* MICA
* MICB
* ULBP1–6

The score was averaged across all cells belonging to the same CAR variant.

## Visualization

A Cell Press-style blue barplot was generated using:

```r
geom_col()
```

## Main outputs

* `CAR_NKG2DL_sum_scores.csv`
* `CAR_total_NKG2DL_score.png`

---

# 4. Correlation between TCAT programs and NKG2D ligands

File:
`Correlation with ligands clean.txt` 

## Purpose

This script evaluates correlations between:

* TCAT T-cell programs
* NKG2D ligand expression

## Method

Spearman correlation analysis was performed between:

* TCAT program usage scores
* individual ligands
* total NKG2D ligand score

Programs related to:

* memory
* exhaustion
* cytotoxicity
* Th subsets
* interferon response
* proliferation
* cell cycle

were analyzed.

## Statistical analysis

Correlations were computed using:

```r
cor.test(method = "spearman")
```

False discovery rate correction was performed using:

```r
p.adjust(method = "BH")
```

## Visualizations

Two major visualizations were generated:

### Heatmap

Correlation coefficient (rho) visualization.

### Dotplot

* color = correlation coefficient
* size = statistical significance

## Main outputs

* `TCAT_vs_NKG2DL_correlation_table_with_SUM.csv`
* `TCAT_vs_NKG2DL_correlation_heatmap_with_SUM.png`
* `TCAT_vs_NKG2DL_correlation_dotplot_with_SUM.png`

---

# 5. Cell cycle and NKG2D ligand analysis

File:
`Cell Cycle.txt` 

## Purpose

This script investigates relationships between:

* cell cycle state
* CAR composition
* NKG2D ligand expression

## Cell cycle scoring

Cell cycle phase scores were calculated using Seurat built-in gene sets:

```r
CellCycleScoring()
```

using:

* S-phase genes
* G2M-phase genes

### S phase

S\ Score

### G2M phase

G2M\ Score

## NKG2D ligand scoring

Per-cell ligand scores were calculated from:

* MICA
* MICB
* ULBP1–6

## Pie-chart visualization

Each CAR was represented as a scatterpie showing:

* TCAT composition
* S phase or G2M score
* NKG2DL score

Coordinates were converted into Z-scores before plotting.

## Main outputs

* `FilteredCAR_PieScatter_S_vs_NKG2DL_with_labels_Zscore.png`
* `FilteredCAR_PieScatter_G2M_vs_NKG2DL_with_labels_Zscore.png`
* `FilteredCAR_PieScatter_S_vs_NKG2DL_NO_LABELS.png`
* `FilteredCAR_PieScatter_G2M_vs_NKG2DL_NO_LABELS.png`

---

# Output summary

The workflow generated:

* filtered Seurat objects
* TCAT annotations
* CAR-specific UMAPs
* TCAT composition analyses
* NKG2D ligand scoring
* cell cycle analyses
* correlation heatmaps
* publication-ready figures

---

# Notes

* Only CAR variants with at least 100 cells were retained for downstream analysis.
* TCAT/starCAT annotations were based on the TCAT.V1 reference.
* All analyses were performed on filtered T-cell populations after exclusion of epithelial/tumor-associated cells.

---

# Citation

If using this workflow, please cite:

> Castellanos-Rueda R, Di Roberto RB, Biebrich F, Schlatter FS et al.
> speedingCARs: accelerating the engineering of CAR T cells by signaling domain shuffling and single-cell sequencing.
> Nat Commun. 2022 Nov 2;13(1):6555.
