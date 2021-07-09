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

# Normalize data
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

# Identify highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# top 10 genes plot
top10 <- head(VariableFeatures(pbmc), 10) # top 10 of the HVG
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
plot2

# Scale data: mu=0, sd=1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear transformation
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) # only for HVGs
# examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # test PCA
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc) # approximating the "true dimensionality / dim rank" of dataset

# -----

pbmc <- FindNeighbors(pbmc, dims = 1:10) # Jaccard similarity: |A,B|/|AvB|
pbmc <- FindClusters(pbmc, resolution = 0.5) # iteratively grouping genes hierarchially
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label=TRUE)

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

# DEG with test.use=ROC
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# Expression overlayed on UMAP
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                               "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Expression heatmap for given cells and features
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assign a priori cell types to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Save R dataset
saveRDS(pbmc, file = "./pbmc3k_final.rds")
sessionInfo()
