Вот собранный **полный README.md** в аккуратном GitHub-формате (можешь просто вставить в репозиторий):

---

# 🧬 TF ChIP-seq Analysis of NKG2D Ligand Promoters

## 📌 Overview

This project performs **transcription factor (TF) binding analysis** at promoter regions of **NKG2D ligand genes** using ENCODE ChIP-seq data.

The goal is to identify TFs that potentially regulate:

* **MICA, MICB**
* **ULBP1–6 (RAET1 family)**

by quantifying TF binding peaks within **promoter regions (TSS ± 2 kb)**.

---

## 📂 Project Structure

```
.
├── data/
│   ├── LIST_of_files.txt        # filtered ENCODE files (.bed.gz only)
│   ├── metadata.tsv             # parsed ENCODE metadata
│   ├── gencode.v49.basic.annotation.gtf
│
├── peaks/
│   ├── *.bed.gz                # TF ChIP-seq peak files
│
├── scripts/
│   ├── 01_metadata_download.sh
│   ├── 02_metadata_parsing.sh
│   ├── 03_promoter_overlap.R
│
├── results/
│   ├── TF_ChIPseq_overlap_long.csv
│   ├── TF_ChIPseq_overlap_matrix.csv
│   ├── heatmaps/
│
└── README.md
```

---

## 🧪 Data Download and Preprocessing

### ENCODE Dataset Selection

TF ChIP-seq datasets were downloaded from the
ENCODE Project

Search filters:

* Organism: **Homo sapiens**
* Assay: **TF ChIP-seq**
* Biosample type: **cell lines**
* Sample type: **cancer cell lines**
* File type: **processed peaks**

### Selected Transcription Factors

```
SP1, TP53, CEBPB, CREB1, SRF, RARA, RXRA, STAT3, RUNX1, SPI1,
TAL1, TCF12, CEBPA, TCF4, POU5F1, TCF3, MYC, RXRB, GATA4,
E2F1, ATF4, FOXO1, SREBF1
```

---

## 📥 Download Procedure

1. ENCODE file list saved to:

   ```
   LIST_of_files.txt
   ```

2. File filtering:

   * Kept only:

     ```
     *.bed.gz
     ```
   * Removed:

     ```
     .bam, .bigWig, .bigBed
     ```

3. File accession extraction from filenames

4. Metadata retrieval via API:

```bash
curl -L https://www.encodeproject.org/files/<FILE_ACCESSION>/?format=json
```

---

## 🧾 Metadata Processing

Metadata fields extracted:

* TF name (`target.label`)
* Cell line (`biosample_ontology.term_name`)
* Assay type
* Lab
* Output type

Final table:

```
metadata.tsv
```

Used to annotate each ChIP-seq peak file.

---

## 🧬 Cell Lines Used

| Cell Line | Number of Experiments |
| --------- | --------------------- |
| K562      | 40                    |
| HepG2     | 33                    |
| A549      | 22                    |
| MCF-7     | 14                    |
| HeLa-S3   | 6                     |
| Ishikawa  | 4                     |
| HCT116    | 3                     |
| SK-N-SH   | 3                     |
| HL-60     | 1                     |
| NB4       | 1                     |

---

## 🧬 Gene Annotation

Annotation used:

```
GENCODE v49 (GRCh38)
```

* File:

  ```
  gencode.v49.basic.annotation.gtf
  ```

### Notes

* Genome build: **hg38**
* Promoter definition:

  ```
  TSS ± 2 kb
  ```
* Alternative GENCODE versions (v38–v44) produce comparable results

---

## 🎯 Target Genes (NKG2D Ligands)

```
MICA, MICB,
ULBP1, ULBP2, ULBP3,
RAET1E (ULBP4),
RAET1G (ULBP5),
RAET1L (ULBP6)
```

---

## 🔬 Analysis Workflow

### Step 1 — Load data

* Peak files (`.bed.gz`)
* Metadata (`metadata.tsv`)
* GTF annotation

### Step 2 — Define promoters

* Extract gene coordinates
* Compute:

```
TSS ± 2 kb
```

### Step 3 — Convert to genomic ranges

Using:

* `GenomicRanges`
* `IRanges`
* `rtracklayer`

### Step 4 — Overlap analysis

* Intersect TF peaks with promoters
* Count overlaps per:

  * TF
  * Gene
  * File

### Step 5 — Aggregation

Generate:

* Long format table
* TF × Gene matrix

---

## 📊 Output Files

### Main results

```
TF_ChIPseq_overlap_long.csv
```

Columns:

* TF
* gene
* cell line
* file accession
* number of peaks

```
TF_ChIPseq_overlap_matrix.csv
```

Matrix:

* rows = TFs
* columns = genes
* values = peak counts

---

## 📈 Visualization

Heatmaps generated using:

* `pheatmap`

Optional normalization:

* peaks per TF
* peaks per file

---

## ⚙️ R Environment

```
R version 4.5.3 (2026-03-11)
Platform: Windows 11 x64
```

### Key packages

```
data.table
dplyr
tidyr
stringr
GenomicRanges
IRanges
rtracklayer
pheatmap
```

---

## 🔁 Reproducibility

To reproduce the analysis:

1. Download ENCODE `.bed.gz` files
2. Generate `metadata.tsv`
3. Download GENCODE v49 annotation
4. Run:

```r
source("scripts/03_promoter_overlap.R")
```

---

## ⚠️ Notes

* Only **processed peak files (.bed.gz)** were used
* Results depend on:

  * ENCODE dataset availability
  * TF coverage per cell line
* Promoter-based analysis is robust to minor annotation differences

---

## 📚 Citation

If using this workflow, cite:

* ENCODE Project
* GENCODE annotation

---

## 👤 Author

**Shynggys Tursymbek**
PhD student — Biomedical Sciences

---

Если хочешь — могу дальше:

* 🔥 добавить **figure section (как в статье Cell/STAR methods)**
* 🔥 сделать **README под Nature / Cell Press стиль**
* 🔥 или сразу собрать **Methods section для статьи из этого README**
