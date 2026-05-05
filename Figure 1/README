## Data Processing and Analysis

The file **"41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"** was obtained from the study:

> Systematic pan-cancer analysis of tumour purity  
> DOI: 10.1038/ncomms9971

Only samples with an **ESTIMATE score ≥ 0.8** were retained and used for downstream analysis.

---

## Data Sources

Gene expression data (**"TcgaTargetGtex_rsem_gene_tpm"**) were downloaded from the UCSC Xena Browser:  
https://xenabrowser.net/

---

## Scripts

### 1. Sample Filtering

**"Script for filtration.R"**  
This script filters the TCGA dataset by removing samples with low tumour purity (ESTIMATE score < 0.8), based on the reference dataset described above.

---

### 2. NKG2D Ligand Analysis

**"TCGA GTEX NKG2D ligands.R"**  
This script performs the analysis of NKG2D ligand expression across TCGA and GTEx samples and generates the final figure used in the article.

Outliers are excluded from visualization in the final plot for clarity; however, they are retained in the dataset and included in all statistical analyses.

---

## Notes

- All analyses were performed in R.
- Raw large-scale datasets (e.g., TCGA/GTEx expression matrices) are not included in this repository due to size limitations.
