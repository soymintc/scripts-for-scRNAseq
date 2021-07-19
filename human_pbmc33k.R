library(dplyr)
library(Seurat)
library(patchwork)

# Load data from single directory
setwd('~/local/human_pbmc33k')
p33k.data = Read10X(data.dir = "./filtered_gene_bc_matrices/hg19")
p33k = CreateSeuratObject(counts = p33k.data, project = "pbmc33k")
p33k[["percent.mt"]] <- PercentageFeatureSet(p33k, pattern = "^MT-")
VlnPlot(p33k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Mean UMI cnt check
umis = p33k@assays$RNA@counts
cell_umi_sum = colSums(umis)
print(sprintf("Mean UMIs per cell: %.2f", mean(colSums(umis))))

# Visualize highly expressed genes
C = thyroid@assays$RNA@counts # counts dataset only
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
        names = C@Dimnames[[1]][most_expressed]) # gene symbols

# Filter out "cells"
pbmc <- subset(pbmc, subset = 
                 nFeature_RNA > 200 & # rm low or high freq expr
                 nFeature_RNA < 2500 & 
                 percent.mt < 5) # rm high mt expr ratio