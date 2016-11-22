
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)

in_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(in_dir,full.names = TRUE)
files <- files[grepl("sort",files) & !grepl("bai",files)]

options(mc.cores = 24)

exo <- list()
exo[[1]] <- ExoData(files[1])
