## Environment

All analyses were performed in the following computational environment:

- **R version:** 4.5.3 (2026-03-11)  
- **Platform:** x86_64-w64-mingw32 (Windows 11 x64)  
- **Time zone:** Asia/Qyzylorda  

---

## Data Sources

### DepMap datasets

The following datasets were used from the Broad DepMap project:

- **OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv**  
  Gene expression matrix (log2(TPM + 1)) for human protein-coding genes across DepMap cell lines.

- **Model.csv**  
  Metadata file containing cell line annotations, including:
  - Cell line names  
  - Tissue lineage  
  - Primary disease classification  

These files are available from the DepMap portal:  
https://depmap.org/portal/download/

---

## Data Processing

- Gene expression values were used directly in log2(TPM + 1) format.  
- Gene names were cleaned to remove annotation suffixes (e.g., `"MICA (XXXX)" → "MICA"`).  
- Duplicate gene columns were removed prior to analysis.  
- Only NKG2D ligand genes were retained:
  - MICA, MICB  
  - ULBP1–3  
  - RAET1E (ULBP4), RAET1G (ULBP5), RAET1L (ULBP6)  

- NKG2D ligand score was calculated as:
  > the mean expression of all detected ligands per cell line  

- Metadata from `Model.csv` was merged using `ModelID`.

---

## Output Files Used in article

The analysis generates:

- `DepMap_NKG2DL_score_distribution_histogram.png`  
  Distribution of NKG2D ligand scores across cell lines  
---

## Key Packages

The following R packages were used across different analysis steps:

### Core data processing
- data.table (v1.18.2.1)  
- dplyr (v1.2.0)  
- stringr (v1.6.0)  

### Data import
- readxl (v1.4.5)  

### Visualization
- ggplot2 (v4.0.2)  

### Annotation (Bioconductor)
- org.Hs.eg.db (v3.21.0)  
- AnnotationDbi (v1.70.0)  

---

## Notes

- Analyses were performed using R-based pipelines for transcriptomic and DepMap data.  
- NKG2D ligand score was calculated as the mean expression of detected ligands per cell line.  
- Visualization was performed using `ggplot2`.  
