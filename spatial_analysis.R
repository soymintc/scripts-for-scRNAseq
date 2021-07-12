library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# install.packages("hdf5r")
# sudo apt-get install libhdf5-dev
library(hdf5r)

setwd('~/local/mouse_sgt_ant1/mouse_brain/')

# Use Load10X_Spatial to load data from the following file structure:
# ├── filtered_feature_bc_matrix
# │   ├── barcodes.tsv.gz
# │   ├── features.tsv.gz
# │   └── matrix.mtx.gz
# ├── filtered_feature_bc_matrix.h5 
# └── spatial 
#     ├── aligned_fiducials.jpg
#     ├── detected_tissue_image.jpg
#     ├── scalefactors_json.json # [req] high -> low-res image scaling info
#     ├── tissue_hires_image.png
#     ├── tissue_lowres_image.png # [req]
#     └── tissue_positions_list.csv # [req] the coordinates for each spot 
brain <- Load10X_Spatial(data.dir="./outs/raw",
                         filename="filtered_feature_bc_matrix.h5",
                         assay='Spatial',
                         filter.matrix = TRUE)

Idents(object = brain) <- "anterior1" # label experiment
Project(object = brain) <- "mouseBrain" # label project

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Filtering cells // change nCount_RNA -> nCount_Spatial (refer Environment)
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")
VlnPlot(brain, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(brain, feature1 = "nCount_Spatial", feature2 = "percent.mt")
plot2 <- FeatureScatter(brain, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1 + plot2 # plot Expr

# Filter out "cells" 
brain <- subset(brain, subset = 
                nFeature_Spatial > 200 & # rm low or high freq expr
                nFeature_Spatial < 8000 & 
                percent.mt < 30) # rm high mt expr ratio

# Normalize data
brain <- NormalizeData(brain, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Identify highly variable genes
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(brain), 10) # top 10 of the HVG
plot1 <- VariableFeaturePlot(brain)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data: mu=0, sd=1
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)

# Visualize expression for selected genes
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1)) # changing transparency
p1 + p2

# Linear transformation
brain <- RunPCA(brain, features = VariableFeatures(object = brain)) # only for HVGs
print(brain[["pca"]], dims = 1:5, nfeatures = 5) # genes per PCs
VizDimLoadings(brain, dims = 1:2, reduction = "pca") # PC amplitude check per gene in PC1-2
DimPlot(brain, reduction = "pca") # PCA for PC1 and PC2
DimHeatmap(brain, dims = 1:15, cells = 500, balanced = TRUE)

# Select significant PCs
brain <- JackStraw(brain, num.replicate = 100)
brain <- ScoreJackStraw(brain, dims = 1:20)
JackStrawPlot(brain, dims = 1:20)
ElbowPlot(brain) # approximating the "true dimensionality / dim rank" of dataset

# brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
# brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
# brain <- FindClusters(brain, verbose = FALSE)
# brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Graph construction & Bottom->up Clustering
brain <- FindNeighbors(brain, dims = 1:20) # Jaccard similarity: |A,B|/|AvB|
brain <- FindClusters(brain, resolution = 0.5) # iteratively grouping genes hierarchially
print(head(Idents(brain), 5))

# Run UMAP
brain <- RunUMAP(brain, dims = 1:20)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# Interactive plot examples
SpatialDimPlot(brain, cells.highlight = 
                 CellsByIdentities(object = brain, 
                                   idents = c(2, 1, 4, 3, 5, 8)), 
               facet.highlight = TRUE, ncol = 3)

SpatialDimPlot(brain, interactive = TRUE)

SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

LinkedDimPlot(brain)

# Find DEGs (only based on inter-cell expression) -> Spatial plot of expressions
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
