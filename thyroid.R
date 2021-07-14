library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd('~/local/seurat_raw_data')
dl <- list.dirs(".", full.names = F, recursive=FALSE)
names(dl) <- c("ATC1","ATC3","PTC1","PTC2","ATC2","PTC3","PTC4","ATC4") # set names to vector
marrow.counts <- Read10X(dl) # names:item vector -> create reads data

thyroid <- CreateSeuratObject(counts = marrow.counts, project = "Thyroid")
thyroid[["percent.mt"]] <- PercentageFeatureSet(thyroid, pattern = "^MT-")

table(thyroid@meta.data$orig.ident) # check overall data

thyroid@meta.data$orig.ident %>%
  head()

VlnPlot(thyroid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

thyroid <- subset(thyroid, subset = 
                 nFeature_RNA > 200 & # rm low or high freq expr
                 nFeature_RNA < 2500 & 
                 percent.mt < 10) # rm high mt expr ratio

# Visualize highly expressed genes
C = thyroid@assays$RNA@counts # counts dataset only
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
        names = C@Dimnames[[1]][most_expressed]) # gene symbols

# What are the ^RPL and ^RPS genes?
rownames(thyroid)[grepl("^RPL", rownames(thyroid))]
rownames(thyroid)[grepl("^RPS", rownames(thyroid))]

# Remove ribosomal protein genes
thyroid = thyroid[!grepl("^RPL", rownames(thyroid)),]
thyroid = thyroid[!grepl("^RPS", rownames(thyroid)),]

# # Normalize data
# thyroid <- NormalizeData(thyroid, 
#                          normalization.method = "LogNormalize",
#                          scale.factor = 10000)
# # Identify highly variable genes
# thyroid <- FindVariableFeatures(thyroid, selection.method = "vst", nfeatures = 2000)
# 
# # top 10 genes plot
# top10 <- head(VariableFeatures(thyroid), 10) # top 10 of the HVG
# plot1 <- VariableFeaturePlot(thyroid, log = TRUE)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2
# 
# # Scale data: mu=0, sd=1
# all.genes <- rownames(thyroid)
# thyroid <- ScaleData(thyroid, features = all.genes)

# # Linear transformation
# thyroid <- RunPCA(thyroid, features = VariableFeatures(object = thyroid)) # only for HVGs
# VizDimLoadings(thyroid, dims = 1:2, reduction = "pca") # test PCA
# DimPlot(thyroid, reduction = "pca")
# DimHeatmap(thyroid, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(thyroid, dims = 1:30, cells = 500, balanced = TRUE)

# SCTransform: replaces "NormalizeData, FindVariableFeatures, ScaleData
thyroid <- SCTransform(thyroid, assay = "RNA", verbose = TRUE, method = "poisson")

# Dim reduction -> Graph construction -> Clustering -> Plotting
thyroid <- RunPCA(thyroid, assay = "RNA", verbose = FALSE)
thyroid <- FindNeighbors(thyroid, reduction = "pca", dims = 1:30)
thyroid <- FindClusters(thyroid, verbose = FALSE)
thyroid <- RunUMAP(thyroid, reduction = "pca", dims = 1:30)

# Plot UMAP by cluster-identity and given-identity
DimPlot(thyroid, reduction = "umap", group.by = c("ident", "orig.ident"))

# # NOTE: This process can take a long time for big datasets, comment out for expediency. More
# # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# # computation time
# thyroid <- JackStraw(thyroid, num.replicate = 100)
# thyroid <- ScoreJackStraw(thyroid, dims = 1:30)
# JackStrawPlot(thyroid, dims = 1:30)
# ElbowPlot(thyroid) # approximating the "true dimensionality / dim rank" of dataset

# # Clustering
# thyroid <- FindNeighbors(thyroid, dims = 1:30) # Jaccard similarity: |A,B|/|AvB|
# thyroid <- FindClusters(thyroid, verbose=TRUE) # iteratively grouping genes hierarchially
# thyroid <- RunUMAP(thyroid, reduction = "pca", dims = 1:30)
# DimPlot(thyroid, reduction = "umap", group.by = c("ident", "orig.ident"))

