library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# install.packages("hdf5r")
# sudo apt-get install libhdf5-dev
library(hdf5r)

setwd('~/local/mouse_sgt_ant1/mouse_brain/outs/filtered')
# Use Load10X_Spatial to load data from the following file structure:
# "ant1" or "post1" directory
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
brain1 <- Load10X_Spatial(data.dir="./ant1",
                         filename="filtered_feature_bc_matrix.h5",
                         assay='Spatial',
                         slice="anterior1",
                         filter.matrix = TRUE)
brain2 <- Load10X_Spatial(data.dir="./post1",
                          filename="filtered_feature_bc_matrix.h5",
                          assay='Spatial',
                          slice="posterior1",
                          filter.matrix = TRUE)

brain1@meta.data$orig.ident = "anterior1" # give cells identity; default: SeuratProject
brain2@meta.data$orig.ident = "posterior1" # give cells identity; default: SeuratProject
Idents(object = brain1) <- "anterior1" # label experiment
Idents(object = brain2) <- "posterior1" # label experiment
Project(object = brain1) <- "brain" # label project
Project(object = brain2) <- "brain" # label project

brain <- merge(brain1, brain2) # summary below:
# An object of class Seurat 
# 32285 features across 6050 samples within 1 assay 
# Active assay: Spatial (32285 features, 0 variable features)

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2) # side-by-side plot

# Filtering cells // change nCount_RNA -> nCount_Spatial (refer Environment)
VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", 
                            "percent_mito", "percent_hb"), 
        pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", 
                                       "percent_mito", "percent_hb"))
# Filter out "cells"
brain = brain[, brain$nFeature_Spatial > 500 & 
                brain$percent_mito < 25 & 
                brain$percent_hb < 20]
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", 
                                       "percent_mito", "percent_hb"))

# Top expressed genes
C = brain@assays$Spatial@counts # counts dataset only
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
        names = C@Dimnames[[1]][most_expressed]) # gene symbols

# Removal of selected genes (e.g. mitochondrial)
dim(brain) # [1] 32285  5806
#   filter Bl1
brain <- brain[!grepl("Bc1", rownames(brain)), ]
#   filter Mitocondrial
brain <- brain[!grepl("^mt-", rownames(brain)), ]
#   filter Hemoglobin gene (optional if Hb genes needed)
brain <- brain[!grepl("^Hb.*-", rownames(brain)), ]
dim(brain) # [1] 32263  5806


# Normalize, Find HVGs, Scale data
brain <- SCTransform(brain, assay = "Spatial", verbose = TRUE, method = "poisson")
#   output as follows:
      # Calculating cell attributes from input UMI matrix: log_umi
      # Variance stabilizing transformation of count matrix of size 19356 by 5806
      # Model formula is y ~ log_umi
      # Get Negative Binomial regression parameters per gene
      # Using 2000 genes, 5000 cells
      # |========================================================================| 100%
      # Found 115 outliers - those will be ignored in fitting/regularization step
      # 
      # Second step: Get residuals using fitted parameters for 19356 genes
      # |========================================================================| 100%
      # Computing corrected count matrix for 19356 genes
      # |========================================================================| 100%
      # Calculating gene attributes
      # Wall clock passed: Time difference of 5.204434 mins
      # Determine variable features
      # Place corrected count matrix in counts slot
      # Centering data matrix
      # |========================================================================| 100%
      # Set default assay to SCT
      # There were 50 or more warnings (use warnings() to see the first 50)

# Marker plot: Hippocampus  & Choroid plexus 
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
SpatialFeaturePlot(brain, features = "Ttr", 
                   pt.size.factor = 1, # point size
                   alpha = c(0.5, 0.8)) # transparency of points

# Dim reduction -> Graph construction -> Clustering -> Plotting
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Plot UMAP, 2 plots grouped by ident and given ident
DimPlot(brain, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain) # cluster plot on low-res image
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(brain), 
               facet.highlight = TRUE, ncol = 5) # 5-col spatial plot

# ----------
# # small test to visualize 2 genes
# SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
# p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
# p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
# p1 + p2
# 
# # Dim reduction: Gene expression vector -> PC vector
# brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
# 
# # Graph construction & Cell(barcode) clustering
# brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
# brain <- FindClusters(brain, verbose = FALSE)
# 
# # Create UMAP cluster for visualization
# brain <- RunUMAP(brain, reduction = "pca", dims = 1:20)
# 
# # cluster 2D plot + spatial plot labeling clusters
# p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
# p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
# p1 + p2




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