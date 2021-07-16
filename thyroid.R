library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Get data in multiple directories
setwd('~/local/seurat_raw_data')
dl <- list.dirs(".", full.names = F, recursive=FALSE)
names(dl) <- c("ATC1","ATC3","PTC1","PTC2","ATC2","PTC3","PTC4","ATC4") # set names to vector
marrow.counts <- Read10X(dl) # names:item vector -> create reads data
thyroid <- CreateSeuratObject(counts = marrow.counts, project = "Thyroid")

thyroid[["percent.mt"]] <- PercentageFeatureSet(thyroid, pattern = "^MT-") # calc read% from mitochondrial genes

table(thyroid@meta.data$orig.ident) # check overall data

VlnPlot(thyroid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# In-filter for cells
thyroid <- subset(thyroid, subset = 
                 nFeature_RNA > 200 & # rm low or high freq expr
                 nFeature_RNA < 2500 & 
                 percent.mt < 25) # rm high mt expr ratio

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

# SCTransform: replaces "NormalizeData, FindVariableFeatures, ScaleData
thyroid <- SCTransform(thyroid, assay = "RNA", verbose = TRUE, method = "poisson")

# Dim reduction -> Graph construction -> Clustering -> Plotting
thyroid <- RunPCA(thyroid, assay = "RNA", verbose = FALSE)
thyroid <- FindNeighbors(thyroid, reduction = "pca", dims = 1:30)
thyroid <- FindClusters(thyroid, verbose = FALSE)
thyroid <- RunUMAP(thyroid, reduction = "pca", dims = 1:30)

# Plot UMAP by cluster-identity and given-identity
DimPlot(thyroid, reduction = "umap", group.by = c("ident", "orig.ident"))
