
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)


dr = "/p/keles/ChIPexo/volume4"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("carro",files)]
files = files[grep("sort",files)]
files = files[grep("mouse",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("txt",files,invert = TRUE)]
files = files[grep("seq",files,invert = TRUE)]

options(mc.cores = 20)

reads = files %>% mclapply(readGAlignments,param = NULL)
reads = reads %>% mclapply(as,"GRanges")
names(reads) = paste("Rep",c(3,1,2))

exo = lapply(reads,function(x)ExoData(reads = x))

as_pbc = exo %>% lapply(function(x)paramDist(x)$beta1)


