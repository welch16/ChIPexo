
rm(list = ls())

library(ChIPexoQual)
library(readr)
library(GenomicAlignments)
library(magrittr)
library(parallel)

dr1 = "/p/keles/ChIPexo/volume4/venters_data/sortbam"
dr2 = "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"

files1 = list.files(dr1,full.names = TRUE)
files1 = files1[grep("TBP",files1)]
files1 = files1[grep("bai",files1,invert = TRUE)]

files2 = list.files(dr2,full.names = TRUE)
files2 = files2[grep("TBP",files2)]
files2 = files2[grep("bai",files2,invert = TRUE)]

files = c(files1,files2)

## blacklist
blacklist = "/p/keles/ChIPexo/volume4/blacklists/wgEncodeHg19ConsensusSignalArtifactRegions.bed"
blacklist = blacklist %>% read_delim(delim = "\t",col_names = FALSE)
blacklist = GRanges(seqnames = blacklist$X1,
                    ranges = IRanges(
                        start = blacklist$X2,
                        end = blacklist$X3))

reads = lapply(files,readGAlignments,param = NULL)
reads = mclapply(reads,as,"GRanges",mc.cores = 5)

## reads in blacklist
ribl = sapply(reads,function(x)x %>% countOverlaps(blacklist) %>%
                               sum)
depth = sapply(reads,length)
## > ribl / depth * 100
## [1] 1.446500 1.023062 1.050551 2.116903 2.126468

exo = lapply(reads,function(x)ExoData(reads = x,mc.cores = 20))

ov = mclapply(exo,findOverlaps,blacklist,mc.cores = 5)

nblacklist = sapply(ov,function(x)x %>% queryHits %>%
                                  unique %>% length)

## > nblacklist
## [1] 18822  4354  3738 25051 20036
## > nblacklist / sapply(exo,length) * 100
## [1] 0.4276918 0.8440242 0.7460393 0.3546692 0.4117999

## type of regions
exo_bl = mapply(function(x,y){
    x[unique(queryHits(y))]
    },exo,ov,SIMPLIFY = FALSE)

## summary for blacklist islads
exo_bl %>% lapply(function(x)x$depth %>% summary)
    ## usually low, it seems there is an outlier region
exo_bl %>% lapply(function(x)x$FSR %>% summary)
    ## in the extremes 0 or 1
exo_bl %>% lapply(function(x)x$uniquePos %>% summary)
    ## usually low too, it seems there is an outlier there
exo_bl %>% lapply(function(x)x$ARC %>% summary)
    ## bad replicates above 1, the other ones below
exo_bl %>% lapply(function(x)x$URC %>% summary)
    ## all close to zero except nexus 1

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
## [1]  28.88414 271.74299 377.96027  18.42574  97.56817
## > beta2_bl
## [1]  0.6689158  8.1245346 10.7900794  1.0214929  5.4548959


library(ggplot2)
library(scales)
figsdr = "figs/NAR_review/blacklist"


plots1 = stats_bl %>% lapply(function(x)
    ggplot(x, aes(uniquePos,depth))+geom_point(alpha = I(1/5))+
    xlim(0,4e3)+ylim(0,6e5))

plots2 = stats_bl %>% lapply(function(x)
    ggplot(x, aes(width,depth))+geom_point(alpha = I(1/5))+
    xlim(0,25e2)+ylim(0,6e5))


pdf(file.path(figsdr,"TBP_scatter_beta1.pdf"))
u = mapply(function(x,z){
    print(x + ggtitle(z))},plots1,c("exo1","exo2","exo3","nexus1","nexus2"))
dev.off()


pdf(file.path(figsdr,"TBP_scatter_beta2.pdf"))
u = mapply(function(x,z){
    print(x + ggtitle(z))},plots2,c("exo1","exo2","exo3","nexus1","nexus2"))
dev.off()


## Exo paramdist without considering those regions

exo_w = mapply(function(x,y){
    x[-queryHits(y)]},exo,ov,SIMPLIFY = FALSE)

## small test
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
## [1]  10.93148 230.70245 274.34997   3.22884  19.69646
## > beta2_w
## [1] -0.03464683  1.80471537  2.17102767  0.02093917  0.13399473

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
## [1]  10.515677 229.195402 268.757758   3.180815  19.526086
## > beta2_rlm
## [1] -0.04044973  1.92699150  2.28882803  0.02126860  0.14151755


## ll = mapply(function(x,y,z){
##     x = x[,repl := y]
##     x = x[,blackl := z]
##     x},c(paramDist_w,paramDist_bl),
##     rep(c("exo1","exo2","exo3","nexus1","nexus2"),2),
##     rep(c("no","yes"),each = 5),SIMPLIFY = FALSE)
## ll = rbindlist(ll)


## pdf(file.path(figsdr,"TBP_boxplots.pdf"))
## ggplot(ll[term == "uniquePos"],
##        aes(blackl,estimate))+geom_boxplot()+
##     facet_grid( . ~ repl)+ylab(expression(beta[1]))+
##     xlab("Overlap with blacklist")+ylim(0,500)
## ggplot(ll[term == "width"],
##        aes(blackl,-estimate))+geom_boxplot()+
##     facet_grid( . ~ repl)+ylab(expression(beta[2]))+
##     xlab("Overlap with blacklist")+ylim(-1,20)
## dev.off()

ll = mapply(function(x,y,z){
    x = x[,repl := y]
    x = x[,blackl := z]
    x},c(paramDist_all,paramDist_w,paramDist_bl),
    rep(c("exo1","exo2","exo3","nexus1","nexus2"),3),
    rep(c("All regions","Don't overlap blacklist","Overlap blacklist"),each = 5),SIMPLIFY = FALSE)
ll = rbindlist(ll)

library(ggplot2)

theme_set(theme_bw())

pdf(file.path(figsdr,"TBP_boxplots.pdf"))
ggplot(ll[term == "uniquePos"],
       aes(blackl,estimate,colour = blackl))+geom_boxplot()+
    facet_grid( . ~ repl)+ylim(0,500)+ylab(expression(beta[1]))+
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

