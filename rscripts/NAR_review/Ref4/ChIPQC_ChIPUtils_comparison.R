
rm(list = ls())

library(ChIPUtils)
library(ChIPQC)
library(parallel)
library(BiocParallel)

dir1 = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
dir2 = "/p/keles/ChIPexo/volume4/venters_data/sortbam"
dir3 = "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"

files1 = list.files(dir1,full.names = TRUE)
files2 = list.files(dir2,full.names = TRUE)
files3 = list.files(dir3,full.names = TRUE)

grepv <- function(pattern,x,invert = FALSE)x[grep(pattern,x,invert = invert)]

files1 = grepv("sort",files1)
files1 = grepv("bai",files1,invert = TRUE)

files2 = grepv("TBP",files2)
files2 = grepv("bai",files2,TRUE)

files3 = grepv("TBP",files3)
files3 = grepv("bai",files3,TRUE)

dr = "/p/keles/ChIPexo/volume3/ChIPexo/otherQC_methods"

foxA1 = list()
foxA1[["FoxA1-rep1"]] = ChIPQCsample(files1[2],annotation = NULL,runCrossCor = TRUE)
foxA1[["FoxA1-rep2"]] = ChIPQCsample(files1[3],annotation = NULL,runCrossCor = TRUE)
foxA1[["FoxA1-rep3"]] = ChIPQCsample(files1[1],annotation = NULL,runCrossCor = TRUE)

save(foxA1,file = file.path(dr,"ChIPQC_FoxA1.RData"))

tbpexo = list()
tbpexo[["TBPexo-rep1"]] = ChIPQCsample(files2[1],annotation = NULL,runCrossCor = TRUE)
tbpexo[["TBPexo-rep2"]] = ChIPQCsample(files2[2],annotation = NULL,runCrossCor = TRUE)
tbpexo[["TBPexo-rep3"]] = ChIPQCsample(files2[3],annotation = NULL,runCrossCor = TRUE)

save(tbpexo,file = file.path(dr,"ChIPQC_TBPexo.RData"))

tbpnexus = list()
tbpnexus[["TBPnexus-rep1"]] = ChIPQCsample(files3[1],annotation = NULL,runCrossCor = TRUE)
tbpnexus[["TBPnexus-rep2"]] = ChIPQCsample(files3[2],annotation = NULL,runCrossCor = TRUE)

save(tbpnexus,file = file.path(dr,"ChIPQC_TBPnexus.RData"))
