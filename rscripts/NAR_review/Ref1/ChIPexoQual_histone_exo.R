
rm( list = ls())

library(ChIPexoQual)
library(parallel)
library(GenomicAlignments)

files = "/p/keles/ChIPexo/volume4/exo_histone_data/BAM/histone1_rep?.sort.bam"
files = Sys.glob(files)

reads = lapply(files,readGAlignments,param = NULL)
