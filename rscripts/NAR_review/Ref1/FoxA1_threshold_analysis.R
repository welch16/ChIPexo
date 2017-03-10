
rm(list = ls())

library(ChIPexoQual)

indir = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(indir,full.names = TRUE)

files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]

library(BiocParallel)
thresh = c(1,5,25,50)

reads = mclapply(files,readGAlignments,param = NULL,mc.cores = 3)

options(mc.cores = 22)

exo = lapply(reads,function(x)ExoData(reads = x ))

island_depth = lapply(exo,function(x)mcols(x)$depth)

## ## respective depth quantiles (calculated with ECDF)
## lapply(island_depth,function(x)ecdf(x)(thresh) * 100)

## [[1]]
## [1]  9.871078 41.904989 98.731624 99.567053

## [[2]]
## [1] 54.36674 95.89683 99.41133 99.68149

## [[3]]
## [1] 60.31944 97.24536 99.71111 99.88926


exoList = list()

exoList[["rep1"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[2]],height = x)})

exoList[["rep2"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[3]],height = x)})

exoList[["rep3"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[1]],height = x)})
    
figsdir = "figs/NAR_review/Threshold_analysis/FoxA1"

## nms = thresh

## exoList = lapply(exoList,
##                  function(x){
##                      names(x) = nms
##                      x})

library(ggplot2)

pdf(paste0(figsdir,"_ARCvURC_threshold.pdf"),height = 5,width = 12)
print(ARCvURCplot(exoList[[1]],names.input = thresh) + ggtitle("Rep-1") + xlim(0,5))
print(ARCvURCplot(exoList[[2]],names.input = thresh) + ggtitle("Rep-2") + xlim(0,5))
print(ARCvURCplot(exoList[[3]],names.input = thresh) + ggtitle("Rep-3") + xlim(0,5))
dev.off()


pdf(paste0(figsdir,"_FSRdist_threshold.pdf"),height = 9,width = 6)
print(FSRDistplot(exoList[[1]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("Rep-1"))
print(FSRDistplot(exoList[[2]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("Rep-2"))
print(FSRDistplot(exoList[[3]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("Rep-3"))
dev.off()


pdf(paste0(figsdir,"_RegionComp_threshold.pdf"),height = 9,width = 6)
print(regionCompplot(exoList[[1]],depth.values = seq_len(50),names.input = thresh) + ggtitle("Rep-1"))
print(regionCompplot(exoList[[2]],depth.values = seq_len(50),names.input = thresh) + ggtitle("Rep-2"))
print(regionCompplot(exoList[[3]],depth.values = seq_len(50),names.input = thresh) + ggtitle("Rep-3"))
dev.off()


tt = as.character(thresh)

pdf(paste0(figsdir,"_beta1_boxplot.pdf"))
print(paramDistBoxplot(exoList[[1]],names.input = tt,sort.as.numeric = TRUE) + ggtitle("Rep-1")+ylim(0,15))
print(paramDistBoxplot(exoList[[2]],names.input = tt,sort.as.numeric = TRUE)+ggtitle("Rep-2")+ylim(0,15))
print(paramDistBoxplot(exoList[[3]],names.input = tt,sort.as.numeric = TRUE)+ggtitle("Rep-3")+ylim(0,15))
dev.off()

pdf(paste0(figsdir,"_beta2_boxplot.pdf"))
print(paramDistBoxplot(exoList[[1]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("Rep-1")+ylim(0,45))
print(paramDistBoxplot(exoList[[2]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("Rep-2")+ylim(0,45))
print(paramDistBoxplot(exoList[[3]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("Rep-3")+ylim(0,45))
dev.off()
