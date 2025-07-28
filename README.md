# Transcriptomic Signatures of the Neocortex in Mesial Temporal Lobe Epilepsy

Mesial Temporal Lobe Epilepsy (mTLE) is a prevalent form of focal epilepsy, traditionally linked to hippocampal dysfunction. However, the neocortex may also contribute significantly to seizure propagation and cognitive impairments. This study leverages RNA sequencing data from the neocortex of mTLE patients and healthy controls to identify differentially expressed genes (DEGs) and explore their associated biological processes. Our analysis revealed 22 significantly upregulated genes enriched in RNA processing and splicing pathways, and 11 downregulated genes associated with synaptic function and vesicle fusion. These findings suggest that neocortical transcriptional changes play a role in mTLE pathology, potentially opening new therapeutic avenues beyond the hippocampus.

This repository contains a complete R-based pipeline for analyzing **differential gene expression** in the **neocortex tissue** of patients with **mesial temporal lobe epilepsy (mTLE)** using RNA-Seq count data from GEO dataset **GSE134697**. The analysis includes preprocessing, DESeq2 statistical testing, visualization (heatmaps and volcano plots), and Gene Ontology (GO) enrichment for biological interpretation.

## Dataset

- **GSE134697_CountMatrix.tsv**: Raw count matrix (downloaded from GEO).
- **GSE134697_metadata.csv**: Sample metadata derived from the GSE134697 series matrix file.
- You can explore the dataset [here on NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134697).

---
# Install base packages
install.packages(c("data.table", "pheatmap", "ggplot2", "dplyr", "gridExtra"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))

---
## Workflow Overview

### 1. **Data Loading & Preparation**
- Load and clean gene expression data.
- Assign sample metadata (condition and tissue type).
- Filter out genes with zero counts across all samples.

### 2. **DESeq2 Differential Expression Analysis**
- Subset to neocortex samples.
- Identify differentially expressed genes (DEGs) using DESeq2.
- Filter DEGs using thresholds: `padj < 0.05` and `|log2FC| > 1`.

### 3. **Visualization**
- **Heatmap**: Top 20 most significant DEGs.
- **Volcano Plot**: Global overview of gene expression changes.

### 4. **GO Enrichment Analysis**
- Identify enriched Biological Process (BP) terms among:
  - Upregulated genes.
  - Downregulated genes.
- Visualized with `dotplot()`.

### 5. **Export Results**
- All DESeq2 results saved as CSV.
- Filtered significant genes saved separately.

---

## R Packages Used

- **[DESeq2]** for differential expression
- **pheatmap** for heatmaps
- **ggplot2** for volcano plots
- **clusterProfiler** for GO analysis
- **gridExtra** for presenting multiple visualizations
- **org.Hs.eg.db** for gene annotations
- **dplyr** and **data.table** for data manipulation

---

## Key Results

- Differential expression analysis between **epileptic vs non-epileptic** neocortex.
- Enriched biological processes hint at molecular pathways involved in mTLE.
- Heatmap and volcano plot clearly show significant gene regulation patterns.

---

## Reproducibility

Ensure the following files are present in your working directory:
- `GSE134697_CountMatrix.tsv`
- `GSE134697_metadata.csv`

### Run the full pipeline:

source("Script.R")

---

## 🔗 Citation

If you use this repository or parts of it, please cite the original dataset and tools used:

> Cagnan H, et al. "GSE134697 - RNA-seq of neocortex and hippocampus from mTLE patients." NCBI GEO (2019).

---

## Author

**Almustafa Zein Elabedein Bashir**
[GitHub: zeinmursal](https://github.com/zeinmursal)
[zein.mursal@gmail.com]

```
