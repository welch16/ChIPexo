
rm( list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4/pugh_data/mace"
files = list.files(dr,full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]

names(files) = paste("Rep",seq_len(3),sep = "-")

library(parallel)

options(mc.cores = 20)

reads = files %>% mclapply(readGAlignments,param = NULL)

reads = reads %>% mclapply(as,"GRanges")

exo = reads %>% lapply(function(x)ExoData(reads = x))
                      

library(dplyr)

build_summary <- function(...,names.input = NULL,summaryFun = mean)
{
    args = list(...) %>% unlist
    if(is.null(names.input))names.input = names(args)

    tibble(
        sample = names.input,
        seqDepth = exo %>% sapply(nreads),
        nIslands = exo %>% sapply(length),
        beta1 = exo %>% sapply(function(x)summaryFun(paramDist(x)$beta1)),
        beta2 = exo %>% sapply(function(x)-summaryFun(paramDist(x)$beta2))
        )
}


## build_summary(exo)
## > build_summary(exo)
## # A tibble: 3 × 5
##   sample seqDepth nIslands    beta1      beta2
##    <chr>    <dbl>    <int>    <dbl>      <dbl>
## 1  Rep-1 23,576,694  5815055 3.141124 0.03504025
## 2  Rep-2 20,947,081  4948524 3.043299 0.03053204
## 3  Rep-3 37,688,587  5840635 5.358850 0.06588083

## > build_summary(exo,summaryFun=median)
## # A tibble: 3 × 5
##   sample seqDepth nIslands    beta1      beta2
##    <chr>    <dbl>    <int>    <dbl>      <dbl>
## 1  Rep-1 23576694  5815055 3.041465 0.03110531
## 2  Rep-2 20947081  4948524 2.963305 0.02695947
## 3  Rep-3 37688587  5840635 5.290756 0.05862988



figs = "figs/NAR_review/CTCF_mace"

pdf(file.path(figs,"CTCF_ARCvURC.pdf"),width = 12)
ARCvURCplot(exo)
ARCvURCplot(exo,both.strand = TRUE)
dev.off()


pdf(file.path(figs,"CTCF_FSRDist.pdf"),width = 6)
FSRDistplot(exo,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(300))
FSRDistplot(exo,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(300),both.strand = TRUE)
dev.off()

pdf(file.path(figs,"CTCF_RegionComp.pdf"),width = 6)
regionCompplot(exo,depth.values = seq_len(50))
dev.off()


pdf(file.path(figs,"CTCF_qc_scores_boxplot.pdf"))
paramDistBoxplot(exo,which.param = "beta1")
paramDistBoxplot(exo,which.param = "beta2")
dev.off()


library(ChIPUtils)

reads_chipu = files %>% mclapply(create_reads)


pbc = reads_chipu %>% sapply(PBC)

## > pbc
##     Rep-1     Rep-2     Rep-3 
## 0.4654122 0.4292791 0.2744155 

library(htSeqTools)


ssd2 = reads %>% sapply(ssdCoverage)
gini2 = reads %>% sapply(giniCoverage,mc.cores = 22)


## > ssd2
##     Rep-1     Rep-2     Rep-3 
## 0.8102357 0.7846769 0.6348150 
## > gini2
##                 Rep-1     Rep-2     Rep-3
## gini        0.9541502 0.9604241 0.9690235
## gini.adjust 0.2183883 0.2029997 0.2513372


library(data.table)

sizes = read.table("/p/keles/SOFTWARE/hg19.chrom.sizes") %>% data.table


scc = lapply(reads_chipu,strand_cross_corr,seq_len(300),sizes,TRUE)


scc = mapply(function(x,y)x[,rep := y],scc,names(exo),SIMPLIFY = FALSE) %>% rbindlist


pdf(file = file.path(figs,"CTCF_SCC.pdf"),height = 5)
ggplot(scc,aes(shift,cross.corr,colour = rep))+geom_line()+
    scale_color_brewer(palette = "Set1")
dev.off()



rl =  reads %>% lapply(width) %>% sapply(median)
## Rep-1 Rep-2 Rep-3 
##    48    48    33 




RSC1 <- function(scc,read_length){
  out <- scc[,max(cross.corr)] / scc[shift == read_length, (cross.corr)]
  return(out)
}

RSC2 <- function(scc,read_length){
  mm <- scc[,min(cross.corr)]
  out <- (scc[,max(cross.corr)] - mm) /( scc[shift == read_length, (cross.corr)] - mm)
  return(out)
}

RSC1(scc[rep == "Rep-1"],rl[1])
RSC1(scc[rep == "Rep-2"],rl[2])
RSC1(scc[rep == "Rep-3"],rl[3])


## > RSC1(scc[rep == "Rep-1"],rl[1])
## [1] 1.260482
## > RSC1(scc[rep == "Rep-2"],rl[2])
## [1] 1.138268
## > RSC1(scc[rep == "Rep-3"],rl[3])
## [1] 1.201708


RSC2(scc[rep == "Rep-1"],rl[1])
RSC2(scc[rep == "Rep-2"],rl[2])
RSC2(scc[rep == "Rep-3"],rl[3])

## > RSC2(scc[rep == "Rep-1"],rl[1])
## [1] 1.275711
## > RSC2(scc[rep == "Rep-2"],rl[2])
## [1] 1.150706
## > RSC2(scc[rep == "Rep-3"],rl[3])
## [1] 1.214258


## > scc[,max(cross.corr)/min(cross.corr),by = rep]
##      rep       V1
## 1: Rep-1 22.82047
## 2: Rep-2 13.79216
## 3: Rep-3 20.51460


