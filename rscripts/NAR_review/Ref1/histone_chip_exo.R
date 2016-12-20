

rm(list = ls())

library(ChIPexoQual)
library(BiocParallel)
library(GenomicAlignments)

indir = "/p/keles/ChIPexo/volume4/exo_histone_data"

files = list.files(indir,full.names = TRUE,recursive = TRUE)
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]


reads = lapply(files,readGAlignments,param = NULL)


files1 = files[grep("histone1",files)]
files2 = files[grep("histone2",files)]
files3 = files[grep("histone3",files)]




exolist1 = lapply(files1,ExoData,mc.cores = 22)
exolist2 = lapply(files2,ExoData,mc.cores = 22)
exolist3 = lapply(files3,ExoData,mc.cores = 22)
