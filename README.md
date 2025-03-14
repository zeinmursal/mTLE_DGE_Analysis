# mTLE_DGE_Analysis

This repository performs a differential gene expression (DGE) analysis on the GSE134697 dataset, focusing on mesial temporal lobe epilepsy (mTLE) samples. The analysis compares gene expression between neocortex and hippocampus samples from 17 patients using R and DESeq2.

## Overview
- **Goal**: Identify differentially expressed genes (DEGs) between neocortex and hippocampus in mTLE patients.
- **Dataset**: [GSE134697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134697) from GEO.
- **Tools**: R, DESeq2, ggplot2, pheatmap, dplyr, tidyr, tibble.
- **Outputs**: Volcano plots, heatmaps, bar plots, and PCA plots for quality control.

## Requirements
- **R**: Version 4.0.0 or higher.
- **R Packages**: 
  ```R
  install.packages(c("readr", "DESeq2", "ggplot2", "pheatmap", "dplyr", "tidyr", "tibble"))
