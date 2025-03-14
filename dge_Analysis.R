# Load required packages
library(readr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)

# Setup Instructions:
# 1. Place this script and GSE134697_CountMatrix.csv in the same directory.
# 2. Set your working directory to this script's location before running:
#    setwd("<path_to_directory_containing_this_script>")
#    Example: setwd("C:/Users/zeinm/Documents/mTLE_DGE_Analysis")
# 3. If the CSV is not included, download it from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134697
#    and update file_path accordingly.

# Load the count matrix with a relative path
file_path <- "GSE134697_CountMatrix.csv"
if (!file.exists(file_path)) {
  stop("Error: GSE134697_CountMatrix.csv not found. Please place the file in the same directory as this script or update the file_path variable with the correct path.")
}
GSE134697_CountMatrix <- read_csv(file_path)

# Prepare count data
all_columns <- colnames(GSE134697_CountMatrix)
sample_names <- all_columns[7:length(all_columns)]  # Exclude first 6 annotation columns
count_data <- as.data.frame(GSE134697_CountMatrix[, 7:ncol(GSE134697_CountMatrix)])
rownames(count_data) <- GSE134697_CountMatrix$Geneid

# Verify data
print(ncol(count_data))  # Should be 36
print(head(count_data))

# Define number of mTLE patients
num_patients_tle <- 17
total_samples <- length(sample_names)

# Create metadata for mTLE samples only
sample_data_tle <- data.frame(
  sample_name = sample_names[1:(2 * num_patients_tle)],
  patient_id = factor(rep(1:num_patients_tle, each = 2)),
  sample_type = factor(rep(c("neocortex", "hippo"), num_patients_tle), levels = c("neocortex", "hippo"))
)
rownames(sample_data_tle) <- sample_data_tle$sample_name

# Subset count data to mTLE samples
count_data_tle <- count_data[, sample_data_tle$sample_name]

# Verify sample names match
if (!all(sample_data_tle$sample_name %in% colnames(count_data_tle))) {
  stop("Sample names in metadata do not match count matrix columns!")
}

# Single DESeq2 analysis with paired design
dds <- DESeqDataSetFromMatrix(countData = count_data_tle, 
                              colData = sample_data_tle, 
                              design = ~ patient_id + sample_type)
dds <- DESeq(dds)

# Extract results for each patient
results_list <- list()
for (k in 1:num_patients_tle) {
  contrast <- c("sample_type", "hippo", "neocortex")
  res_k <- results(dds, contrast = contrast, 
                   lfcThreshold = 0,  # No fold change threshold here, handle in Significant
                   alpha = 0.1)  # Relaxed to 0.1 for more DEGs
  
  # Add patient-specific identifier
  res_k$patient <- k
  results_list[[k]] <- res_k
  
  # Volcano Plot (relaxed thresholds)
  res_df <- as.data.frame(res_k) %>%
    mutate(Gene = rownames(res_k),
           Significant = padj < 0.1 & abs(log2FoldChange) > 0.5) %>%  # Relaxed thresholds
    replace_na(list(Significant = FALSE))
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste("Volcano Plot - Patient", k),
         x = "Log2 Fold Change (Hippo vs Neocortex)",
         y = "-log10 Adjusted P-value") +
    theme_minimal()
  ggsave(paste0("output/volcano_patient_", k, ".png"), p, width = 8, height = 6)
}

# Genes significant in 5+ patients (relaxed thresholds)
significant_genes <- lapply(results_list, function(res) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>%  # Relaxed thresholds
    pull(gene)
})
gene_counts <- table(unlist(significant_genes))
genes_of_interest <- names(gene_counts)[gene_counts >= 5]
print(genes_of_interest)

# Heatmap
if (length(genes_of_interest) > 0) {
  count_subset <- count_data_tle[genes_of_interest, sample_data_tle$sample_name]
  dds_all <- DESeqDataSetFromMatrix(countData = count_subset, 
                                    colData = sample_data_tle, 
                                    design = ~ sample_type)
  vst <- assay(vst(dds_all))
  annotation_col <- data.frame(SampleType = sample_data_tle$sample_type, 
                               Patient = sample_data_tle$patient_id,
                               row.names = sample_data_tle$sample_name)
  pheatmap(vst, 
           annotation_col = annotation_col,
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           show_rownames = length(genes_of_interest) <= 50,
           main = "Heatmap of Genes Significant in 5+ Patients",
           filename = "output/heatmap_significant_genes.png",
           width = 10, 
           height = 8)
}

# Bar Plot for Documentation (relaxed thresholds)
deg_counts <- sapply(results_list, function(res) {
  sum(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5, na.rm = TRUE)  # Relaxed thresholds
})
deg_df <- data.frame(
  Patient = factor(1:num_patients_tle),
  DEG_Count = deg_counts
)

bar_plot <- ggplot(deg_df, aes(x = Patient, y = DEG_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Number of Significant DEGs per Patient",
       x = "Patient ID",
       y = "Number of DEGs (padj < 0.1, |log2FC| > 0.5)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/deg_counts_per_patient.png", bar_plot, width = 10, height = 6)

# Check Data Quality, Generate a PCA Plot
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "sample_type")
ggsave("output/pca_plot.png")