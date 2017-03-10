
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)
library(dplyr)
library(readr)

dr = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(dr,full.names = TRUE,pattern = "sort.bam")
files = files[grep("bai",files,invert = TRUE)]

names(files) = paste("Rep",c(3,1,2),sep = "-")

reads = files %>% mclapply(readGAlignments,param = NULL) %>%
    mclapply(as,"GRanges")

exo = reads %>% lapply(function(x)ExoData(reads = x))

## read blacklist

blacklist = "/p/keles/ChIPexo/volume4/blacklists/mm9-blacklist.bed"
blacklist = blacklist %>% read_delim(delim = "\t",col_names = FALSE)
blacklist = GRanges(seqnames = blacklist$X1,
                    ranges = IRanges(
                        start = blacklist$X2,
                        end = blacklist$X3))


ribl = sapply(reads,function(x)x %>% countOverlaps(blacklist) %>%
                               sum)
depth = sapply(reads,length)
## ribl / depth * 100
## [1] 2.472307 1.818175 2.077000

ov = mclapply(exo,findOverlaps,blacklist,mc.cores = 3)

nblacklist = sapply(ov,function(x)x %>% queryHits %>%
                                  unique %>% length)
## > nblacklist
## [1]  6059 10988 12923
## > nblacklist / sapply(exo,length) * 100
## [1] 0.2151421 0.1296940 0.1101174

## type of regions

exo_bl = mapply(function(x,y){
    x[unique(queryHits(y))]
    },exo,ov,SIMPLIFY = FALSE)

## summary for blacklist islads
exo_bl %>% lapply(function(x)x$depth %>% summary)
    ## usually low, except one outlier
exo_bl %>% lapply(function(x)x$FSR %>% summary)
    ## in the extremes 0 or 1
exo_bl %>% lapply(function(x)x$uniquePos %>% summary)
    ## usually low too, except one outlier
exo_bl %>% lapply(function(x)x$ARC %>% summary)
    ## below .5, except the outlier
exo_bl %>% lapply(function(x)x$URC %>% summary)
    ## close to one except rep3


source("~/Desktop/Docs/Code/ChIPexoQual/R/base_summaryStats.R")
source("~/Desktop/Docs/Code/ChIPexoQual/R/base_summaryPlots.R")

ntimes = 1000
nregions = 1000

## small test
library(data.table)
stats_bl = exo_bl %>% lapply(function(x)x %>% as.data.frame %>%
                                  as.data.table)
stats_bl = stats_bl %>% lapply(function(x)x[,.(depth,width,uniquePos)])

paramDist_bl = stats_bl %>%
    lapply(function(x){
        pp = mclapply(seq_len(ntimes),
                      calculateParamDist,x,nregions,mc.cores =10)
        do.call(rbind,pp)})

beta1_bl = paramDist_bl %>%
    sapply(function(x)x[term == "uniquePos",median(estimate)])
beta2_bl = paramDist_bl %>%
    sapply(function(x)x[term == "width",-median(estimate)])

## > beta1_bl
## [1] 62.411032  7.438226  6.575949
## > beta2_bl
## [1] 4.7781657 0.7759577 0.8017310

## library(ggplot2)
## library(scales)
## figsdr = "figs/NAR_review/blacklist"


## plots1 = stats_bl %>% lapply(function(x)
##     ggplot(x, aes(uniquePos,depth))+geom_point(alpha = I(1/5))+
##       xlim(0,1e3)+ylim(0,25e3))

## plots2 = stats_bl %>% lapply(function(x)
##     ggplot(x, aes(width,depth))+geom_point(alpha = I(1/5))+
##     ylim(0,25e3))

## pdf(file.path(figsdr,"FoxA1_scatter_beta1.pdf"))
## u = mapply(function(x,z){
##     print(x + ggtitle(z))},plots1,c("rep3","rep1","rep2"))
## dev.off()


## pdf(file.path(figsdr,"FoxA1_scatter_beta2.pdf"))
## u = mapply(function(x,z){
##     print(x + ggtitle(z))},plots2,c("rep3","rep1","rep2"))
## dev.off()



## Exo paramdist without considering those regions

exo_w = mapply(function(x,y){
    x[-queryHits(y)]},exo,ov,SIMPLIFY = FALSE)

## without blacklist
library(data.table)
stats_w = exo_w %>% lapply(function(x)x %>% as.data.frame %>%
                                  as.data.table)
stats_w = stats_w %>% lapply(function(x)x[,.(depth,width,uniquePos)])

paramDist_w = stats_w %>%
    lapply(function(x){
        pp = mclapply(seq_len(ntimes),
                      calculateParamDist,x,nregions,mc.cores =10)
        do.call(rbind,pp)})

beta1_w = paramDist_w %>%
    sapply(function(x)x[term == "uniquePos",median(estimate)])
beta2_w = paramDist_w %>%
    sapply(function(x)x[term == "width",-median(estimate)])


## > beta1_w
## [1] 8.125377 1.934753 1.505312
## > beta2_w
## [1] 0.046611535 0.017712744 0.009387129


## all

stats_all = exo %>% lapply(function(x)x %>% as.data.frame %>%
                                      as.data.table)

stats_all = stats_all %>% lapply(function(x)x[,.(depth,width,uniquePos)])

paramDist_all = stats_all %>%
    lapply(function(x){
        pp = mclapply(seq_len(ntimes),
                      calculateParamDist,x,nregions,mc.cores =10)
        do.call(rbind,pp)})

beta1_all = paramDist_all %>%
    sapply(function(x)x[term == "uniquePos",median(estimate)])
beta2_all = paramDist_all %>%
    sapply(function(x)x[term == "width",-median(estimate)])


blacklist_scores  = mapply(function(x,y,z){
    x = x[,repl := y]
    x = x[,blackl := z]
    x},c(paramDist_all,paramDist_w,paramDist_bl),
    rep(c("Rep-3","Rep-1","Rep-2"),3),
    rep(c("All regions","Don't overlap blacklist","Overlap blacklist"),each = 3),SIMPLIFY = FALSE) %>%
    rbindlist

FSR = exo %>% lapply(.FSRDistDataFrame,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(300),FALSE)

RegionComp = exo %>% lapply(.regionCompDataFrame,depth.values = seq_len(50))

ARCvURC = exo %>% lapply(function(x){   
    x %>% mcols %>% as.data.frame %>% as.tbl %>%
        mutate(seqnames = as.character(seqnames(x)) ,
               start = start(x),
               end = end(x),
               width = width(x)) %>%
        select(seqnames,start,end,width,everything())})

combine_tbl <- function(x){
    mapply(function(z,y){z %>% mutate(Repl = y)},x,names(x),SIMPLIFY = FALSE) %>%
        bind_rows}

FSR = FSR %>% lapply(function(x){x %>% as.data.frame %>% as.tbl})
RegionComp = RegionComp %>% lapply(function(x){x %>% as.data.frame %>% as.tbl})

ll = ll %>% as.tbl

ARCvURC = combine_tbl(ARCvURC)
FSR = combine_tbl(FSR)
RegionComp = combine_tbl(RegionComp)


outdr = "data/figures/fig3"

write_tsv(ll, path = file.path(outdr,"fig3D_blacklist.tsv"))

write_tsv(ARCvURC, path = file.path(outdr,"fig3A_ARCvURC.tsv"))


write_tsv(FSR, path = file.path(outdr,"fig3B_FSR.tsv"))
write_tsv(RegionComp, path = file.path(outdr,"fig3C_RegionComp.tsv"))
