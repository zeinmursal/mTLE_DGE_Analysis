# mTLE_DGE_Analysis
Differential gene expression (DGE) analysis on the GSE134697 dataset, comparing neocortex and hippocampus samples in mesial temporal lobe epilepsy (mTLE) using R and DESeq2.

## Overview
- **Goal**: Identify DEGs in mTLE patients.
- **Dataset**: [GSE134697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134697).
- **Tools**: R, DESeq2, ggplot2, pheatmap.

## Setup
1. Clone: `git clone https://github.com/zeinmursal/mTLE_DGE_Analysis.git`
2. Place `GSE134697_CountMatrix.csv` in `data/` (download from GEO).
3. Install R packages: `install.packages(c("readr", "DESeq2", "ggplot2", "pheatmap", "dplyr", "tidyr", "tibble"))`
4. Set working directory: `setwd("C:/path/to/mTLE_DGE_Analysis")`
5. Knit: `rmarkdown::render("scripts/mTLE_DGE_Analysis.Rmd")`

## Results
See [`mTLE_DGE_Analysis.html`](mTLE_DGE_Analysis.html) for the full report. Sample outputs in `output/`:
- `volcano_patient_1.png`
- `heatmap_significant_genes.png`
- `deg_counts_per_patient.png`
- `pca_plot.png`

## License
[MIT License](LICENSE)
