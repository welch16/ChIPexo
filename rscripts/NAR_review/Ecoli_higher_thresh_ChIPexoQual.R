
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(parallel)

dr = "/p/keles/ChIPexo/volume4/landick_data/ChIPexo"
files = list.files(dr,full.names = TRUE,recursive = TRUE)
files = files[grep("txt",files,invert = TRUE)]
files = files[grep("sort",files)]
files = files[grep("bam",files)]
files = files[grep("bai",files,invert =TRUE)]

sig70 = as.character(c(931,933,935,937,1311,1314,1317,1320))

files = files[sapply(sig70,grep,files)]

thresh = c(1,3,5,10,25,50,75,100)
cores = 10

reads = mclapply(files,readGAlignments,param = NULL,mc.cores = cores)
reads = mclapply(reads,as,"GRanges",mc.cores = cores)

ntimes = 1e3

exo = lapply(thresh,
             function(z){
                 lapply(reads,function(x,tt)ExoData(reads  =x ,height = tt,mc.cores= cores,ntimes = ntimes),z)
             })

names(exo) = as.character(thresh)

library(magrittr)

nms = files %>% basename %>% strsplit("_") %>% sapply(function(x)x[1])

exo = lapply(exo,function(x){
    names(x) = nms
    x})

exo = exo %>% unlist

exo_param = lapply(exo,paramDist)

library(dplyr)
library(tidyr)

exo_param = mapply(function(x,y){
    x %>% as.data.frame %>% as.tbl %>% mutate(nm = y)
    },exo_param,names(exo_param),SIMPLIFY = FALSE)

exo_param = do.call(rbind,exo_param) %>% separate(nm,into = c("thresh","edsn"),sep ="\\.",)

tt = as.character(thresh)

exo_param = exo_param %>%
    mutate(thresh = factor(thresh,levels = tt))

figsdr = "figs/NAR_review/Threshold_analysis/Ecoli"

library(ggplot2)

pdf(file.path(figsdr,"Sig70_scores.pdf"),width = 12,height = 5)
ggplot(exo_param,aes(thresh,beta1))+geom_boxplot()+
    facet_grid( . ~ edsn)+theme(axis.text.x = element_text(angle = 90))+ylab(expression(beta[1]))
ggplot(exo_param,aes(thresh,-beta2))+geom_boxplot()+
    facet_grid( . ~ edsn)+theme(axis.text.x = element_text(angle = 90))+ylab(expression(beta[2]))
dev.off()

idx = exo %>% names %>% strsplit("\\.")
idx = do.call(rbind,idx) %>% as.data.frame %>% as.tbl
names(idx) = c("thresh","edsn")


pdf(file.path(figsdr,"Sig70_ARCvURC.pdf"),width = 8)
lapply(nms,
       function(x){
           ARCvURCplot(exo[idx$edsn ==x],names.input = thresh)+ggtitle(x)+xlim(0,4)
       })
dev.off()


pdf(file.path(figsdr,"Sig70_FSRDist.pdf"),height = 10)
lapply(nms,
       function(x){
           FSRDistplot(exo[idx$edsn == x],names.input = thresh,depth.values = seq_len(100),
                       quantiles = c(.1,.25,.5,.75,.9))+ggtitle(x)})
dev.off()

pdf(file.path(figsdr,"Sig70_RegionComp.pdf"),height = 10)
lapply(nms,
       function(x){
           regionCompplot(exo[idx$edsn == x],names.input = thresh,depth.values = seq_len(30))+ggtitle(x)})
dev.off()
