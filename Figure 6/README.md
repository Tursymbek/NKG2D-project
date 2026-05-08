# Taiji-based transcription factor network analysis

## Overview

To identify transcriptional regulatory programs associated with high NKG2D ligand (NKG2DL) expression, a multi-omics network analysis workflow was performed using the Taiji pipeline (v1.3.0). The analysis integrated chromatin accessibility (ATAC-seq peak regions) with RNA expression data from cancer cell lines exhibiting elevated NKG2DL expression signatures.

The workflow aimed to:

* identify transcription factors (TFs) potentially regulating NKG2DL-associated programs,
* construct TF-centered regulatory networks,
* detect highly connected hub TFs,
* identify transcriptional modules using graph-based clustering approaches,
* and prioritize key regulatory nodes associated with immune signaling, stress responses, differentiation, and nuclear receptor pathways.

---

# Software and computational environment

## Taiji

* Version: Taiji v1.3.0
* Website: Taiji pipeline ([web-site](https://taiji-pipeline.github.io/), [git-hub](https://github.com/Taiji-pipeline/Taiji))
* Executed in:

  * Ubuntu (WSL2)
  * Windows 11 host system

## Additional R packages

The downstream network analysis and visualization were performed in R using:

* `igraph`
* `ggraph`
* `ggplot2`
* `ggforce`
* `data.table`

---

# Genome annotation and reference files

## GENCODE annotation

The following genome annotation file was used:

```text id="t9x0uk"
gencode.v49.basic.annotation.gtf
```

This annotation file corresponds to:

* Human genome assembly: hg38 / GRCh38
* Source: [GENCODE Release 49](https://www.gencodegenes.org/human/release_49.html), [download gencode.v49.basic.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz)

The GTF annotation was used by Taiji to:

* define gene coordinates,
* assign ATAC-seq peaks to nearby genes,
* generate promoter-associated regulatory regions,
* and connect transcription factor binding events to target genes.

---

# Input files used for Taiji

The project structure contained the following directories:

```text id="6jstzv"
Taiji_Project/
├── config/
├── RNA/
├── ATAC/
├── meta/
└── results/
```

---

# RNA expression input

RNA expression quantification tables were stored in:

```text id="z6h8v8"
RNA/
```

The following files were used:

```text id="j0wkfi"
A549_gene_quant.tsv
HEPG2_gene_quant.tsv
NCIH929_gene_quant.tsv
PANC1_gene_quant.tsv
PC3_gene_quant.tsv
SKNSH_gene_quant.tsv
```

These files contained normalized gene-level expression values for each cell line and were used by Taiji to:

* infer transcriptional activity,
* estimate TF regulatory influence,
* and integrate expression with chromatin accessibility data.

---

# ATAC-seq input

ATAC peak files were stored in:

```text id="d6sg7g"
ATAC/
```

The following peak files were used:

```text id="j1wd4v"
A549_peak_rep1.bed.gz
HEPG2_peak_rep1.bed.gz
NCIH929_peak_rep1.bed.gz
PANC1_peak_rep1.bed.gz
PC3_peak_rep1.bed.gz
SKNSH_peak_rep1.bed.gz
```

These BED files contained accessible chromatin regions identified by ATAC-seq.

The files were used to:

* identify open chromatin regions,
* predict transcription factor binding sites,
* generate TF–target regulatory edges,
* and construct regulatory networks.

---

# Sample metadata

The metadata table:

```text id="7gk5vx"
meta/samples.tsv
```

was used to define:

* sample identifiers,
* matching RNA and ATAC files,
* and experimental grouping.

---

# Taiji configuration

The workflow was configured using:

```text id="2pg5gc"
config/config.yaml
```

The configuration file defined:

* genome assembly (`hg38`),
* annotation file paths,
* sample names,
* ATAC-seq peak files,
* RNA expression files,
* output directories,
* and multithreading settings.

---

# Taiji workflow execution

The pipeline was executed using:

```bash
/home/tursymbek/taiji run --config config/config.yaml -n 6 +RTS -N6
```

The workflow performed:

* genome indexing,
* ATAC peak processing,
* motif scanning,
* TF binding site prediction,
* TF-target linkage generation,
* expression integration,
* and regulatory network construction.

---

# Important Taiji processing steps

The workflow included several major computational stages:

| Step                   | Purpose                                |
| ---------------------- | -------------------------------------- |
| `ATAC_Merge_Peaks`     | Merge accessible chromatin peaks       |
| `ATAC_Find_TFBS_Union` | Predict TF binding motifs within peaks |
| `Create_Linkage`       | Link TF binding sites to target genes  |
| `RNA_Make_Expr_Table`  | Build expression matrices              |
| `Output_Ranks_SC`      | Generate TF ranking outputs            |

---

# STRING network analysis

To further investigate regulatory relationships among transcription factors, TFs identified from Taiji analysis were imported into the STRING database for protein–protein interaction (PPI) analysis.

STRING database:
[STRING database](https://string-db.org)

---

# STRING input gene list

The TF network contained 105 connected transcription factors. 

These genes were used for:

* network construction,
* hub identification,
* modularity analysis,
* and MCL clustering.

---

# STRING files used

The following STRING export files were used:

```text id="d3pk4o"
string_interactions.tsv
string_node_degrees.tsv
string_MCL_clusters.tsv
```

---

# Purpose of each STRING file

| File                      | Purpose                           |
| ------------------------- | --------------------------------- |
| `string_interactions.tsv` | Protein–protein interaction edges |
| `string_node_degrees.tsv` | Degree centrality information     |
| `string_MCL_clusters.tsv` | MCL cluster assignments           |

---

# Network filtering

The STRING interaction network was imported into R using `igraph`.

Disconnected components and isolated nodes were removed by retaining only the largest connected component of the network.

This filtering step reduced noise and allowed downstream analyses to focus on the core transcriptional regulatory architecture.

---

# Hub transcription factor analysis

Network centrality metrics were calculated, including:

* degree centrality,
* betweenness centrality,
* closeness centrality,
* and PageRank.

An integrated centrality score was generated from normalized centrality metrics to prioritize key transcriptional regulators.

The top-ranked hub transcription factors included:

* TP53,
* STAT3,
* SP1,
* RUNX1,
* RXRA,
* and CEBPB.

These TFs occupied central positions within the regulatory network and connected multiple functional transcriptional modules.

---

# MCL clustering analysis

The transcription factor interaction network was partitioned using the Markov Cluster Algorithm (MCL), a graph-based clustering approach widely used for protein interaction networks.

MCL clustering was preferred over k-means clustering because:

* it directly operates on graph topology,
* detects densely connected modules,
* and is more appropriate for network-structured biological data.

The resulting network modules corresponded to distinct biological programs, including:

* nuclear receptor signaling,
* stress-response transcription,
* immune-associated regulation,
* differentiation-associated TF programs,
* and cell cycle transition modules.

---

# Visualization

Network visualization was performed in R using:

* `ggraph`
* `ggplot2`
* `ggforce`

The final visualizations included:

* MCL cluster coloring,
* hub TF highlighting,
* node scaling based on network degree,
* and cloud-style cluster annotations.

Top hub TFs were highlighted in red, while cluster regions were visualized as semi-transparent modular clouds.
