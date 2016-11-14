
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)
library(dplyr)

## This code shows:
##
## 1 - The pipeline is really easy to run
## 2 - There is tiny efficiency gain when using the reads parameter.
## 3 - 

options(mc.cores = 24)

file <- "/p/keles/ChIPexo/volume4/meijsing_data/K562_GR_chip-exo.sort.bam"
reads <- readGAlignments(file,param = NULL)

sample.depth <- seq(20,50,by = 10) * 1e6
exoList <- ExoDataSubsampling(reads = reads,sample.depth = sample.depth,
                              nregions = 1e3,ntimes = 1e3,
                              verbose = TRUE)

qc_scores <- lapply(exoList,paramDist)

save(qc_scores,file = "K562_GR_parameters_QC.RData")

qc_scores <- mapply(function(x,y){
    x$depth = y
    x },qc_scores,sample.depth,SIMPLIFY = FALSE)


QC <- do.call(rbind,qc_scores) %>% as.data.frame %>% as.tbl

figs_dir <- "figs/NAR_review"

library(ggplot2)
library(scales)

summar <- QC %>% group_by(depth) %>%
    summarise(beta1 = mean(beta1),
              beta2 = -mean(beta2))

theme_set(theme_bw())

pdf(file = file.path(figs_dir,"GR_K562_qc_scores.pdf"))
summar %>%
    ggplot(aes(depth,beta1))+geom_line()+geom_point(size = 3)+ylim(0,10)
summar %>%
    ggplot(aes(depth,beta2))+geom_line()+geom_point(size = 3)+ylim(0,.5)
dev.off()



