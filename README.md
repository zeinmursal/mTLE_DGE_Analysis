# Enhanced Differential Gene Expression Analysis of Mesial Temporal Lobe Epilepsy in Neocortex Using Synthetic Controls

This repository contains an R-based pipeline for analyzing differential gene expression (DGE) in neocortex tissue of patients with mesial temporal lobe epilepsy (mTLE) using RNA-Seq count data from GEO dataset GSE134697. The study explores the impact of limited control sample sizes by incorporating synthetic control samples generated via negative binomial modeling. The pipeline includes preprocessing, DESeq2 statistical testing, visualization (PCA, volcano plots, heatmaps, and GO dotplots), and Gene Ontology (GO) enrichment for biological interpretation.

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
- Load and clean gene expression data from GSE134697 (17 mTLE, 2 control neocortex samples).
- Assign sample metadata (condition and tissue type).
- Filter out genes with total raw counts < 10 across all samples.

### 2. **Synthetic Control Generation**
- Generate three synthetic control samples using negative binomial distribution, with mean and dispersion estimated from the two real controls using DESeq2.
- Example R function for synthetic control generation:
  ```{r}
    simulate_synthetic_control <- function(control_counts, dispersions, sample_name) {
    synthetic_control <- numeric(nrow(control_counts))
    names(synthetic_control) <- rownames(control_counts)
    for (i in seq_len(nrow(control_counts))) {
      mu <- mean(control_counts[i, ], na.rm = TRUE)
      if (is.na(mu) || mu <= 0) mu <- 1
      disp <- ifelse(i <= length(dispersions), dispersions[i], 0.1)
      size_param <- 1/disp
      synthetic_control[i] <- rnbinom(1, mu = mu, size = size_param)
    }
    return(synthetic_control)
  }
  ```

### 3. **DESeq2 Differential Expression Analysis**
- Perform DGE analysis using DESeq2 for two scenarios:
  - Real-only: 17 mTLE vs. 2 controls (identified 33 DEGs: 26 up-regulated, 7 down-regulated).
  - Synthetic-augmented: 17 mTLE vs. 5 controls (2 real + 3 synthetic, identified 1863 DEGs: 1510 up-regulated, 353 down-regulated).
- Filter DEGs using thresholds: padj < 0.05 and |log2FC| > 1.

### 4. **Visualization**
- PCA Plot: Shows clearer separation between mTLE and control samples in synthetic-augmented analysis.
- Volcano Plot: Displays global gene expression changes.
- Heatmap: Top variable genes in synthetic-augmented analysis.
- GO Dotplot: Visualizes enriched biological processes.
- Log2FC Distribution: Boxplot and histogram of log2 fold changes, highlighting improved variance modeling and sensitivity in synthetic-augmented analysis.

### 5. **GO Enrichment Analysis**
- Use clusterProfiler to identify enriched Biological Process (BP) terms for:
- Upregulated DEGs (e.g., synaptic signaling: GO:0007268, with genes like SNAP25, SYT1, GRIN2B).
- Downregulated DEGs.
- Top enriched GO terms (synthetic-augmented):
- Synaptic signaling (GO:0007268, p=1.2e-05, 45 genes).
- Inflammatory response (GO:0006954, p=3.4e-04, 30 genes).
- Signal transduction (GO:0007165, p=5.6e-04, 25 genes).
- Immune response (GO:0006955, p=8.9e-04, 20 genes).
- Intracellular signaling (GO:0035556, p=1.1e-03, 18 genes).

### 5. **Export Results**
- Save DESeq2 results and filtered significant DEGs as CSV files.
- Visualizations saved as PNG files (e.g., pca_plot_synthetic.png, volcano_plot_synthetic.png).
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

- Synthetic augmentation significantly increased DEG detection (33 to 1863), highlighting the impact of sample size on statistical power.
- Enriched pathways include synaptic signaling, neuroinflammation, and immune activation, with upregulated genes like GFAP, C1QA, and SERPINA3 aligning with known mTLE pathology.
- Synthetic controls improved variance modeling, as shown by a narrower log2FC interquartile range and a median log2FC of 1.742, indicating sensitivity to moderate expression changes.
- Note: Results are preliminary due to synthetic controls; biological interpretations require caution and experimental validation.

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

> Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.

> Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012;16(5):284-287.

---

## Author

**Almustafa Zein Elabedein Bashir**
[GitHub: zeinmursal](https://github.com/zeinmursal)
[zein.mursal@gmail.com]

## About
This repository supports an exploratory RNA-seq study focusing on the use of synthetic controls to enhance DGE analysis in underpowered datasets, using mTLE as a case study. It is intended for educational and methodological purposes, not definitive biological claims.
```
