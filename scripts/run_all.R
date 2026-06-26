# -----------------------------------------------------------------------------
# NKG2D-project analysis script order
# -----------------------------------------------------------------------------
# These scripts reproduce the generated source-data files and figures when the
# corresponding public raw datasets are available locally. Several steps require
# large external files (TCGA/GTEx Xena matrix, DepMap 26Q1, ENCODE peak files,
# GSE214231 10x matrices, GSE249511 Salmon quant files) and should be run inside
# their matching analysis folder or with paths configured via environment variables.
#
# Run scripts manually in this order:
# 01_tcga_gtex_nkg2dl_score.R
# 02_depmap_nkg2dl_expression_distribution.R
# 03_depmap_nkg2dl_crispr_gene_effect.R
# 04_taiji_tf_activity_summary.R
# 05_encode_tf_chipseq_promoter_overlap.R
# 06_gse214231_starcat_data_preparation.R
# 07_gse214231_umap_and_composition.R
# 08_gse214231_nkg2dl_car_scores.R
# 09_gse214231_cell_cycle_nkg2dl.R
# 10_gse214231_marker_dotplot.R
# 11_gse214231_tcat_nkg2dl_correlation.R
# 12_gse214231_ligand_correlations_clean.R
# 13_gse249511_bulk_rnaseq_gene_annotation.R
# 14_gse249511_bulk_dotplot_with_nkg2dl.R
# -----------------------------------------------------------------------------
