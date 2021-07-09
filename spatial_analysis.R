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

# Linear transformation
brain <- RunPCA(brain, features = VariableFeatures(object = brain)) # only for HVGs
print(brain[["pca"]], dims = 1:5, nfeatures = 5) # genes per PCs
VizDimLoadings(brain, dims = 1:2, reduction = "pca") # PC amplitude check per gene in PC1-2
DimPlot(brain, reduction = "pca") # PCA for PC1 and PC2
DimHeatmap(brain, dims = 1:15, cells = 500, balanced = TRUE)



# Visualize expression for selected genes
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
