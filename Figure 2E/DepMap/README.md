# NKG2D Ligand Expression-Based Stratification of Cancer Cell Lines

## Overview

This project quantifies **NKG2D ligand (NKG2DL) expression** across selected cancer cell lines using DepMap RNA-seq data and stratifies them into **HIGH** and **LOW** groups based on total ligand expression for following Taiji (v1.3) analysis.

The analysis is fully reproducible and implemented in R.

---

## Biological Rationale

NKG2D ligands (MICA, MICB, and ULBP family members) are key determinants of immune recognition by cytotoxic lymphocytes. Their expression varies substantially across cancer cell lines and reflects underlying cellular stress, immune signaling, and metabolic states.

To capture this variability, ligand expression is aggregated into a single quantitative score per cell line.

---

## NKG2D Ligand Score

The total NKG2D ligand expression score is calculated as:

"NKG2DL_{sum} = \sum_{i=1}^{n} Expression_{ligand_i}"

### Ligands included:

* MICA
* MICB
* ULBP1
* ULBP2
* ULBP3
* RAET1E (ULBP4)
* RAET1G (ULBP5)
* RAET1L (ULBP6)

---

## Data Sources

* **DepMap RNA-seq dataset**
  `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`

* **DepMap metadata**
  `Model.csv`

---

## ENCODE Data Filtering Criteria

ATAC-seq datasets were obtained from ENCODE using strict filtering criteria to ensure consistency and biological relevance.

### Selection filters:

* **Assay type:** ATAC-seq
* **Organism:** *Homo sapiens*
* **Biosample classification:** cell line
* **Cell type:** cancer cell
* **Biosample treatment:** no treatment

---

## Rationale for Filtering

These filters were applied to:

* ensure all samples represent **untreated baseline chromatin states**
* avoid confounding effects from:

  * drug treatments
  * cytokine stimulation
  * genetic perturbations
* maintain **comparability across cell lines**

---

## Additional Constraints

Following filtering, only a limited number of cancer cell lines were available in ENCODE:

* HCT116
* K562
* A549
* HepG2
* MCF-7
* NCI-H929
* PC-3
* Panc1

This limitation directly influenced downstream cell line selection.

---

## Analysis Workflow

### 1. Cell line mapping

Cell line names were matched using `Model.csv`:

* `StrippedCellLineName` used for consistent mapping
* Corresponding `ModelID` extracted

### 2. RNA-seq filtering

* Expression matrix filtered by `ModelID`
* Only **default RNA-seq entries** retained
* Duplicate entries removed

### 3. Data transformation

* Expression matrix transposed to:

  * rows = genes
  * columns = selected cell lines
* Gene names cleaned (removed Ensembl annotations)

### 4. Ligand extraction

Subset of genes corresponding to NKG2D ligands extracted.

### 5. Score calculation

For each cell line:

* expression values summed across ligands
* resulting in **NKG2DL_sum**

### 6. Stratification

Cell lines ranked by NKG2DL_sum and grouped into:

* **HIGH expression**
* **LOW expression**

---

From this set, a subset of biologically relevant cell lines was selected for downstream analysis and stratification based on NKG2D ligand expression:

| Cell Line | Group |
| --------- | ----- |
| A549      | High  |
| HCT116    | High  |
| PC-3      | High  |
| HepG2     | Low   |
| NCI-H929  | Low   |

---

## Code

Main processing script:



---

## Software Environment

R session used:



### Key packages:

* data.table (1.18.2.1)
* dplyr (1.2.0)
* stringr (1.6.0)

---

## Output Files

* `DepMap_selected_cell_lines_RNAseq.csv`
  Full expression matrix (gene × cell line)

* `DepMap_selected_cell_lines_mapping.csv`
  Mapping between cell line names and ModelID

* `DepMap_selected_cell_lines_NKG2DL.csv`
  Expression of NKG2D ligands only

---

## Important Notes

* DepMap expression values are provided as **log2(TPM + offset)**
* Analysis is performed directly on these values for relative comparison
* NKG2DL score reflects **relative ligand abundance**, not absolute expression

---

## Limitations

* Binary grouping (HIGH vs LOW) simplifies continuous variation

* Differences between cell lines may reflect:

  * tissue origin
  * genetic background
  * global transcriptional programs

* The analysis identifies **associations**, not causality

---

## Summary

This workflow provides a reproducible approach to:

* quantify NKG2D ligand expression across cancer cell lines
* compare ligand abundance
* stratify models for downstream immunological analysis

---

## Author

Shynggys Tursymbek
PhD Candidate, Biomedical Sciences
