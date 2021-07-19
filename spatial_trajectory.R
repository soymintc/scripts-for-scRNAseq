# load packages
library(SPATA)
library(magrittr)
library(ggplot2)

# ggplot2 plots can be easily combined with 'patchwork'
library(patchwork)

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
brain1 = initiateSpataObject_10X(input_path="./ant1",
                                 sample_names = "brain1")

# plot a cluster feature
p1 <-
  plotSurface(object = brain1,
              of_sample = "",
              color_to = "seurat_clusters",
              pt_clrp = "npg",
              pt_size = 1.5) + labs(color = "Clusters")

# plot gene-set expression 
p2 <-
  plotSurface(object = brain1,
              of_sample = "",
              color_to = "Ttr",
              pt_size = 1.5,
              smooth = TRUE,
              pt_clrsp = "magma")

# combine with patchwork 
p1 + theme(legend.position = "top") +
  p2 + theme(legend.position = "top") # relocate the legend as the gene set name is long

# interactive plots
plots <- plotSurfaceInteractive(brain1) 

# Surface plots for various genes
plotSurfaceComparison(object = brain1,
                      variables = c("Ttr", "Ppp1r1b", "Gpr88", "Penk", "Pde10a", "Nrgn"),
                      smooth = TRUE,
                      smooth_span = 0.02,
                      pt_size = 1,
                      pt_clrsp = "inferno")

# Expression across subgroups
# e.g. violinplot 
genes_of_interest = c("Ttr", "Hpca", "Nrgn", "Gpr88")
plotDistributionAcross(object = brain1,
                       variables = genes_of_interest[1:3],
                       across = "seurat_clusters",
                       plot_type = "violin",
                       clrp = "npg",
                       ncol = 1) +
  theme(legend.position = "none")


# for understanding purpose

plotSurface(object = brain1,
            color_to = "seurat_clusters",
            pt_clrp = "npg") +
  theme(legend.position = "top") +
  labs(color = "Seurat Clusters")


# Distrubition of gene expr across subtypes: ridgeplot
plotDistributionAcross(object = brain1,
                       variables = genes_of_interest[1:4],
                       across = "seurat_clusters",
                       plot_type = "ridgeplot",
                       nrow = 2) +
  theme(legend.position = "none")


# Manual segmentation 
createSegmentation(brain1)

# display the current segmentation
plotSegmentation(brain1, pt_size = 2.5)
