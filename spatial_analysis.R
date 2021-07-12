library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# install.packages("hdf5r")
# sudo apt-get install libhdf5-dev
library(hdf5r)

setwd('~/local/mouse_sgt_ant1/mouse_brain/')
dl <- list.dirs(".", full.names = F, recursive=FALSE)

# Use Load10X_Spatial to load data from the following file structure:
# ├── filtered_feature_bc_matrix
# │   ├── barcodes.tsv.gz
# │   ├── features.tsv.gz
# │   └── matrix.mtx.gz
# ├── filtered_feature_bc_matrix.h5 
# └── spatial 
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

# # Filtering cells // change nCount_RNA -> nCount_Spatial (refer Environment)
# brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")
# VlnPlot(brain, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) 
# plot1 <- FeatureScatter(brain, feature1 = "nCount_Spatial", feature2 = "percent.mt")
# plot2 <- FeatureScatter(brain, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
# plot1 + plot2 # plot Expr

# # Filter out "cells" 
# brain <- subset(brain, subset = 
#                 nFeature_Spatial > 200 & # rm low or high freq expr
#                 nFeature_Spatial < 8000 & 
#                 percent.mt < 30) # rm high mt expr ratio

# Normalize, Find HVGs, Scale data
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# small test to visualize 2 genes
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

# Dim reduction: Gene expression vector -> PC vector
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

# Graph construction & Cell(barcode) clustering
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)

# Create UMAP cluster for visualization
brain <- RunUMAP(brain, reduction = "pca", dims = 1:20)

# cluster 2D plot + spatial plot labeling clusters
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2




# de_markers <- FindMarkers(brain, ident.1 = 2, ident.2 = 3)
# SpatialFeaturePlot(object = brain, 
#                    features = rownames(de_markers)[1:3], 
#                    alpha = c(0.1, 1), ncol = 3)
# 
# brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", 
#                                        features = VariableFeatures(brain)[1:1000],
#                                        selection.method = "markvariogram")
# 
# brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", 
#                                        selection.method = "markvariogram")
# 
# 
# top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 
#                      6)
# SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))