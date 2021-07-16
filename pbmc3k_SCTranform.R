library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./local/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter out "cells"
pbmc <- subset(pbmc, subset = 
                 nFeature_RNA > 200 & # rm low or high freq expr
                 nFeature_RNA < 2500 & 
                 percent.mt < 5) # rm high mt expr ratio

# SCTransform instead of Seurat normalization of PBMC
pbmc = SCTransform(pbmc, assay="RNA", verbose=TRUE, method="poisson")

# Linear transformation
pbmc <- RunPCA(pbmc, assay="SCT") # using SCTransform matrix

# Draw graph -> Cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10) # Jaccard similarity: |A,B|/|AvB|
pbmc <- FindClusters(pbmc, resolution = 0.5) # iteratively grouping genes hierarchially

# UMAP dimension reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
new.colors = c("#f8766d", "#d39200", "#93aa00", 
               "#00ba38", "#555555", "#619cff", "#00b9e3", "#00c19f", "#db72fb", "#ff61c3")
new.labelcolors = c(rep("black", 0), "red", rep("black", 9))
DimPlot(pbmc, reduction = "umap", label=TRUE,
        cols=new.colors, label.color=new.labelcolors) # one cluster unknown

# Finding DEGs (cluster biomarkers)
# find all markers of cluster 2
# min.pct 
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 1, 2, 3, 4), min.pct = 0.25)
head(cluster5.markers, n = 50)

# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
cluster_markers = pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
cluster_markers

# DEG with test.use=ROC
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("GZMK", "CCL5"), slot="counts", log=T)
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE) # raw counts

# Expression overlayed on UMAP
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                               "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Expression heatmap for given cells and features
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene, group.colors = new.colors) + NoLegend()

# Assign a priori cell types to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "UNKNOWN", "NK", "FCGR3A+ Mono", "CD8 T", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids) # replace index to name vector value

# UMAP after label replacement
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, 
        cols = new.colors, label.color = new.labelcolors) + NoLegend()

# Find genes differently expressed in UNKNOWN compared to CD8 and NK
cluster4.markers = FindMarkers(pbmc, ident.1 = "UNKNOWN", ident.2 = c("CD8 T", "NK"), min.pct = 0.25)
head(cluster4.markers)
VlnPlot(pbmc, features = c("NKG7", "GZMB", "LTB"))
