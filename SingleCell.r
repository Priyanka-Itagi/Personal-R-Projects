# Install and loading Seurat
if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")
library(Seurat)
library(ggplot2)
library(dplyr)

# single-cell RNA-seq data (example dataset)
data <- Read10X(data.dir = "C:\\Users\\Lenovo\\Desktop\\single.csv")
seurat_obj <- CreateSeuratObject(counts = data)

# Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data
seurat_obj <- ScaleData(seurat_obj)

# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Cluster cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Visualize clusters
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1.5) +
  ggtitle("UMAP Visualization of Cell Clusters") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
print(p)

# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save results
write.csv(markers, "single_cell_markers.csv")

# Visualize top markers for a specific cluster (e.g., cluster 0)
top_markers <- markers %>%
  filter(cluster == 0) %>%
  top_n(10, avg_log2FC)

p <- FeaturePlot(seurat_obj, features = top_markers$gene, pt.size = 1.5) +
  ggtitle("Top Marker Genes for Cluster 0") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
print(p)