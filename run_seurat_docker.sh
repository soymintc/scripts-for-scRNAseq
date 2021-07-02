docker run -it \
   -e ROOT=true \
   -e PASSWORD=Tjdals12#$ \
   -v $PWD:/home/rstudio/materials \
   -p 8789:8787 \
   --name sm_seurat \
   satijalab/seurat
