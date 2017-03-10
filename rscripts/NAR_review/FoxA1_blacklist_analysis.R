
rm(list = ls())

library(ChIPexoQual)
library(readr)
library(GenomicAlignments)
library(magrittr)
library(parallel)

dr = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(dr,full.names = TRUE)
files = files[grep("bam",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]

blacklist = "/p/keles/ChIPexo/volume4/blacklists/mm9-blacklist.bed"
blacklist = blacklist %>% read_delim(delim = "\t",col_names = FALSE)
blacklist = GRanges(seqnames = blacklist$X1,
                    ranges = IRanges(
                        start = blacklist$X2,
                        end = blacklist$X3))

reads = mclapply(files,readGAlignments,param = NULL,mc.cores = 3)
reads = mclapply(reads,as,"GRanges",mc.cores = 3)

## reads in blacklist
ribl = sapply(reads,function(x)x %>% countOverlaps(blacklist) %>%
                               sum)
depth = sapply(reads,length)
## ribl / depth * 100
## [1] 2.472307 1.818175 2.077000

exo = lapply(reads,function(x)ExoData(reads = x,mc.cores = 20))

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

library(ggplot2)
library(scales)
figsdr = "figs/NAR_review/blacklist"


plots1 = stats_bl %>% lapply(function(x)
    ggplot(x, aes(uniquePos,depth))+geom_point(alpha = I(1/5))+
      xlim(0,1e3)+ylim(0,25e3))

plots2 = stats_bl %>% lapply(function(x)
    ggplot(x, aes(width,depth))+geom_point(alpha = I(1/5))+
    ylim(0,25e3))

pdf(file.path(figsdr,"FoxA1_scatter_beta1.pdf"))
u = mapply(function(x,z){
    print(x + ggtitle(z))},plots1,c("rep3","rep1","rep2"))
dev.off()


pdf(file.path(figsdr,"FoxA1_scatter_beta2.pdf"))
u = mapply(function(x,z){
    print(x + ggtitle(z))},plots2,c("rep3","rep1","rep2"))
dev.off()



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

## > beta1_all
## [1] 8.082664 1.928376 1.525872
## > beta2_all
## [1] 0.04471741 0.01811462 0.01025854

## calculateParamDist <-
## function(i,stats,nregions)
## {
##     dt <- stats[sample(.N,nregions)]
##     model <- MASS:::rlm(depth ~ 0 + uniquePos + width , data = dt)
##     data.table(broom::tidy(model))
## }


## paramDist_rlm = stats_w %>%
##     lapply(function(x){
##         pp = mclapply(seq_len(ntimes),
##                       calculateParamDist,x,nregions,mc.cores =10)
##         do.call(rbind,pp)})

## beta1_rlm = paramDist_rlm %>%
##     sapply(function(x)x[term == "uniquePos",median(estimate)])
## beta2_rlm = paramDist_rlm %>%
##     sapply(function(x)x[term == "width",-median(estimate)])

## > beta1_rlm
## [1] 8.139362 1.958852 1.513907
## > beta2_rlm
## [1] 0.045892553 0.018607094 0.009722554

## boxplots

ll = mapply(function(x,y,z){
    x = x[,repl := y]
    x = x[,blackl := z]
    x},c(paramDist_all,paramDist_w,paramDist_bl),
    rep(c("rep-3","rep-1","rep-2"),3),
    rep(c("All regions","Don't overlap blacklist","Overlap blacklist"),each = 3),SIMPLIFY = FALSE)
ll = rbindlist(ll)

library(ggplot2)

theme_set(theme_bw())

pdf(file.path(figsdr,"FoxA1_boxplots.pdf"))
ggplot(ll[term == "uniquePos"],
       aes(blackl,estimate,colour = blackl))+geom_boxplot()+
    facet_grid( . ~ repl)+ylim(0,100)+ylab(expression(beta[1]))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")
ggplot(ll[term == "uniquePos"],
       aes(blackl,estimate,colour = blackl))+geom_boxplot()+
    facet_grid( . ~ repl)+ylab(expression(beta[1]))+scale_y_log10()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")
ggplot(ll[term == "width"],
       aes(blackl,-estimate,colour = blackl))+geom_boxplot()+
    facet_grid( . ~ repl)+ylab(expression(beta[2]))+ylim(-.1,5)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")
dev.off()

