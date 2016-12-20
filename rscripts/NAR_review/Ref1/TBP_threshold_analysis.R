
rm(list = ls())

library(ChIPexoQual)

indir = "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files = list.files(indir,full.names = TRUE)

files = files[grep("TBP",files)]
files = files[grep("bai",files,invert = TRUE)]

library(BiocParallel)
thresh = c(1,5,25,50)

reads = lapply(files,readGAlignments,param = NULL)

exoList = list()

exoList[["rep1"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[1]],height = x,mc.cores = 4)})

exoList[["rep2"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[2]],height = x,mc.cores = 4)})

exoList[["rep3"]] = lapply(thresh,function(x){
    message("Using threshold: ",x)
    ExoData(reads = reads[[3]],height = x,mc.cores = 4)})
    
figsdir = "figs/NAR_review/Threshold_analysis/TBPexo"

## nms = thresh

## exoList = lapply(exoList,
##                  function(x){
##                      names(x) = nms
##                      x})


library(ggplot2)

pdf(paste0(figsdir,"_ARCvURC_threshold.pdf"),height = 5,width = 12)
print(ARCvURCplot(exoList[[1]],names.input = thresh) + ggtitle("rep1") + xlim(0,5))
print(ARCvURCplot(exoList[[2]],names.input = thresh) + ggtitle("rep2") + xlim(0,5))
print(ARCvURCplot(exoList[[3]],names.input = thresh) + ggtitle("rep3") + xlim(0,5))
dev.off()


pdf(paste0(figsdir,"_FSRdist_threshold.pdf"),height = 9,width = 6)
print(FSRDistplot(exoList[[1]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("rep1"))
print(FSRDistplot(exoList[[2]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("rep2"))
print(FSRDistplot(exoList[[3]],depth.values = seq_len(120),quantiles = c(.1,.25,.5,.75,.9),names.input = thresh) + ggtitle("rep3"))
dev.off()


pdf(paste0(figsdir,"_RegionComp_threshold.pdf"),height = 9,width = 6)
print(regionCompplot(exoList[[1]],depth.values = seq_len(50),names.input = thresh) + ggtitle("rep1"))
print(regionCompplot(exoList[[2]],depth.values = seq_len(50),names.input = thresh) + ggtitle("rep2"))
print(regionCompplot(exoList[[3]],depth.values = seq_len(50),names.input = thresh) + ggtitle("rep3"))
dev.off()

tt = as.character(thresh)

pdf(paste0(figsdir,"_beta1_boxplot.pdf"))
print(paramDistBoxplot(exoList[[1]],names.input = tt,sort.as.numeric = TRUE) + ggtitle("rep1"))
print(paramDistBoxplot(exoList[[2]],names.input = tt,sort.as.numeric = TRUE)+ggtitle("rep2"))
print(paramDistBoxplot(exoList[[3]],names.input = tt,sort.as.numeric = TRUE)+ggtitle("rep3"))
dev.off()

pdf(paste0(figsdir,"_beta2_boxplot.pdf"))
print(paramDistBoxplot(exoList[[1]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("rep1"))
print(paramDistBoxplot(exoList[[2]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("rep2"))
print(paramDistBoxplot(exoList[[3]],which.param = "beta2",names.input = tt,sort.as.numeric = TRUE)+ggtitle("rep3"))
dev.off()
