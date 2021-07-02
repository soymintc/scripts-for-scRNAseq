docker run -p 8788:8787 \
    -e ROOT=TRUE \
    -e PASSWORD=Tjdals12#$ \
    -v $(pwd):/home/rstudio/local \
    -v /zdisk/kyeonghunjeong/:/zdisk/kyeonghunjeong/ \
    --name sm_rstudio4 \
    rocker/rstudio 
