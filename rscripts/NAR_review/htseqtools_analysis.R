
rm(list = ls())

library(htSeqTools)
library(GenomicAlignments)
library(parallel)
library(BiocParallel)
library(magrittr)

## FoxA1 in mouse liver
indir = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files1 = list.files(indir,full.names = TRUE)
files1 = files1[grep("sort",files1)]
files1 = files1[grep("bai",files1,invert = TRUE)]

reads = lapply(files1,readGAlignments, param = NULL)
reads = reads %>% lapply(as,"GRanges")

ssd1 = reads %>% sapply(ssdCoverage)
gini1 = reads %>% sapply(giniCoverage,mc.cores = 22)

rm(reads)

## ChIP-exo: TBP in human K562
indir = "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files2 = list.files(indir,full.names = TRUE)
files2 = files2[grep("TBP",files2)]
files2 = files2[grep("bai",files2,invert = TRUE)]

reads = lapply(files2,readGAlignments,param = NULL)
reads = reads %>% lapply(as,"GRanges")

ssd2 = reads %>% sapply(ssdCoverage)
gini2 = reads %>% sapply(giniCoverage,mc.cores = 22)

rm(reads)

## ChIP-nexus: TBP in human K562
indir = "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
files3 = list.files(indir,full.names = TRUE)
files3 = files3[grep("K562",files3)]
files3 = files3[grep("bai",files3,invert = TRUE)]

reads = lapply(files3,readGAlignments,param = NULL)
reads = reads %>% lapply(as,"GRanges")

ssd3 = reads %>% sapply(ssdCoverage)
gini3 = reads %>% sapply(giniCoverage,mc.cores = 22)

rm(reads)

htseqtools = list(foxA1 = list(files = files1,ssd = ssd1 , gini = gini1),
           tbp_exo = list(files = files2,ssd = ssd2 , gini = gini2),
           tbp_nexus = list(files = files3,ssd = ssd3 , gini = gini3))

save(htseqtools , file = "otherQC_methods/htseqtools/demo_run.RData")


