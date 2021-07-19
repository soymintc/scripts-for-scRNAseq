# sudo apt-get install libfftw3-dev libfftw3-doc # bash
# sudo apt-get install r-base-dev libcurl4-openssl-dev libxml2-dev libssl-dev libmagick++-6.q16-dev
# sudo apt install -y libproj-dev

base::install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')

install.packages("sf")
install.packages("fftwtools")
BiocManager::install("EBImage")

devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATA")

library(SPATA)
