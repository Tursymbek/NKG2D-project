## DepMap CRISPR Gene Effect Analysis

### Data Sources

The following datasets were used from the Broad DepMap project:

- **CRISPRGeneEffect.csv**  
  Genome-wide CRISPR-Cas9 gene knockout screening data (Chronos gene effect scores), where:
  - 0 → non-essential gene  
  - -0.5 → moderate dependency  
  - -1 → strong dependency (pan-essential range)  

- **Model.csv**  
  Metadata file containing cell line annotations.

Data are available from:  
https://depmap.org/portal/download/

---

## Data Processing

- Gene names were cleaned to remove annotation suffixes (e.g., `"MICA (XXXX)" → "MICA"`).  
- Duplicate gene columns were removed prior to analysis.  
- Only NKG2D ligand genes were analyzed:
  - MICA, MICB  
  - ULBP1–3  
  - RAET1E (ULBP4), RAET1G (ULBP5), RAET1L (ULBP6)  

- DepMap gene naming was harmonized for visualization:
  - RAET1E → ULBP4  
  - RAET1G → ULBP5  
  - RAET1L → ULBP6  

- Metadata from `Model.csv` was merged using `ModelID`.

---

## Analysis

- Gene effect scores were analyzed across all DepMap cell lines.  
- Distributions were evaluated using:
  - Boxplots (per gene)  
  - Density plots  
  - Global histograms  

- Dependency thresholds were indicated:
  - 0 (no effect)  
  - -0.5 (selective dependency)  
  - -1 (common essential genes)  

- Summary statistics were calculated per gene:
  - Mean and median gene effect  
  - Range (min–max)  
  - Percentage of cell lines below dependency thresholds  

---

## Output Files that used in Article

The analysis generates:

- `DepMap_NKG2DL_CRISPRGeneEffect_density.png`  
  Density distribution of gene effects  

- `DepMap_NKG2DL_CRISPRGeneEffect_histogram.png`  
  Overall distribution of gene effects  
