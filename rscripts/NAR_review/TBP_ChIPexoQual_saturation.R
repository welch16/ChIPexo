rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("TBP",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("chipseq",files,invert = TRUE)]

library(parallel)


set.seed(12345)

samp =c( seq_len(4) * 5e6 , seq(3,5) * 1e7)
lvs = c("5,000,000","10,000,000","15,000,000","20,000,000",
                                          "30,000,000","40,000,000","50,000,000")

## files = files[-4]
## reads = lapply(files,readGAlignments,param = NULL)
## reads = lapply(reads,as,"GRanges")



## demo
reads = readGAlignments(files[4],param = NULL)
reads = as(reads,"GRanges")

base = ExoData(reads = reads,mc.cores = 24)


idd = sample(c(TRUE,FALSE),length(base),replace = TRUE)

sample(length(base) / 2,length(base))

base0 = base[!idd]
base1 = base[idd]

groups = 10


idd_reads = sample(seq_len(groups),length(reads),replace = TRUE)


samp_reads = lapply(seq_len(groups),
                    function(i){                        
                        idx = which(idd_reads <= i)
                        reads[idx]})
                       
## saturation
exo = lapply(samp_reads,function(r)ExoData(reads = r , mc.cores = 24))

ov = mclapply(exo,findOverlaps,base,mc.cores = 22)


prob = sapply(ov,function(x){
    a = subjectHits(x) %>% unique %>% length
    a / subjectLength(x)})

library(dplyr)
library(ggplot2)

figs = "figs/NAR_review/threshold"

dt = tibble(reads = samp_reads %>% sapply(length),
            prob)

pdf(file.path(figs,"demo_ChIPnexus1_saturation.pdf"))
dt %>% ggplot(aes(reads,prob))+geom_line()+ylim(0,1)+
    scale_x_continuous(breaks = c(0,1e7,2e7,3e7),
                       labels = c("0","10,000,000","20,000,000","30,000,000"))+
    geom_vline(xintercept = 2e7,linetype = 2,colour = "red")
dev.off()






## demo
reads = readGAlignments(files[5],param = NULL)
reads = as(reads,"GRanges")

base = ExoData(reads = reads,mc.cores = 24)


groups = 40


idd_reads = sample(seq_len(groups),length(reads),replace = TRUE)


samp_reads = lapply(seq_len(10),
                    function(i){                        
                        idx = which(idd_reads <= i)
                        reads[idx]})
                       
## saturation
exo = lapply(samp_reads,function(r)ExoData(reads = r , mc.cores = 24))

ov = mclapply(exo,findOverlaps,base,mc.cores = 22)


prob = sapply(ov,function(x){
    a = subjectHits(x) %>% unique %>% length
    a / subjectLength(x)})

library(dplyr)
library(ggplot2)

figs = "figs/NAR_review/threshold"

dt = tibble(reads = samp_reads %>% sapply(length),
            prob)

pdf(file.path(figs,"demo_ChIPnexus2_saturation.pdf"))
dt %>% ggplot(aes(reads,prob))+geom_line()+ylim(0,1)+
    scale_x_continuous(breaks = c(0,1e7,2e7,3e7),
                       labels = c("0","10,000,000","20,000,000","30,000,000"))+
    geom_vline(xintercept = 2e7,linetype = 2,colour = "red")
dev.off()
