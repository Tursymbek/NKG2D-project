# NKG2D-project

Analysis repository for **Toward fratricide-resistant NKG2D receptor-based CAR-T therapy**.

This version reorganizes the deposited GitHub archive so that folder names match the manuscript STAR Methods and KEY RESOURCES TABLE.

## Repository structure

```text
scripts/                                      Custom R scripts, renamed in analysis order
analyses/                                     Generated outputs by STAR Methods section
  01_TCGA_GTEx_transcriptomic_analysis/
  02_DepMap_expression_and_gene_dependency/
  03_Taiji_transcription_factor_network/
  04_ENCODE_ChIPseq_promoter_overlap/
  05_GSE214231_single_cell_starCAT/
  06_GSE249511_bulk_RNAseq/
Supplementary_Materials/                      Source-data package matching KRT/Table S1-S8 labels
Supplementary_Table_Map.csv                   Manuscript-to-file map
environment/                                  R session information recovered from archive
data/external/README.md                       Instructions for large public raw data
```

## Script order

Run them manually after downloading the required raw public datasets:

1. `scripts/01_tcga_gtex_nkg2dl_score.R`
2. `scripts/02_depmap_nkg2dl_expression_distribution.R`
3. `scripts/03_depmap_nkg2dl_crispr_gene_effect.R`
4. `scripts/04_taiji_tf_activity_summary.R`
5. `scripts/05_encode_tf_chipseq_promoter_overlap.R`
6. `scripts/06_gse214231_starcat_data_preparation.R`
7. `scripts/07_gse214231_umap_and_composition.R`
8. `scripts/08_gse214231_nkg2dl_car_scores.R`
9. `scripts/09_gse214231_cell_cycle_nkg2dl.R`
10. `scripts/10_gse214231_marker_dotplot.R`
11. `scripts/11_gse214231_tcat_nkg2dl_correlation.R`
12. `scripts/12_gse214231_ligand_correlations_clean.R`
13. `scripts/13_gse249511_bulk_rnaseq_gene_annotation.R`
14. `scripts/14_gse249511_bulk_dotplot_with_nkg2dl.R`

`script_manifest.csv` records the original archive filename for each script.

## Analyses matched to STAR Methods

| STAR Methods section | Repository folder |
|---|---|
| TCGA and GTEx transcriptomic analysis of NKG2D ligand expression | `analyses/01_TCGA_GTEx_transcriptomic_analysis/` |
| DepMap analysis of NKG2D ligand expression and gene dependency | `analyses/02_DepMap_expression_and_gene_dependency/` |
| Transcription factor regulatory network analysis using Taiji | `analyses/03_Taiji_transcription_factor_network/` |
| ENCODE ChIP-seq integration with NKG2D ligand promoter regions | `analyses/04_ENCODE_ChIPseq_promoter_overlap/` |
| Single-cell RNA-seq preprocessing and starCAT annotation / TCAT correlations / CAR composition and cell-cycle analyses | `analyses/05_GSE214231_single_cell_starCAT/` |
| Bulk RNA-seq reanalysis of NKG2D-based CAR T cells | `analyses/06_GSE249511_bulk_RNAseq/` |

## External data

Large raw public datasets are intentionally not committed. See `data/external/README.md` for download targets and expected filenames.
