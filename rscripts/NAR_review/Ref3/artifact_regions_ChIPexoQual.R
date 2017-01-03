
rm(list = ls())

library(magrittr)
library(dplyr)
library(readr)
library(data.table)
library(ChIPexoQual)

dr = "/p/keles/ChIPexo/volume4/meijsing_data"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

xlsfiles = files[grep("xls",files)]
readfiles = files[grep("bam",files)]


## load("data/ChIPexo_QC_runs/meijsing_GR_K562_Rep1.RData")

peaks = lapply(xlsfiles,read_delim,delim = "\t",comment = "# ")
peaks[[1]] = peaks[[1]][-1,]
peaks[[3]] = peaks[[3]][-1,]

names(peaks) = c("Input","Negative","Peaks")

library(GenomicRanges)

## stats = ext_stats$stats

## statsgr = stats[,GRanges(seqnames = seqnames,
##                        ranges = IRanges(
##                            start = start,
##                            end = end))]

peaksgr = lapply(peaks,function(x){
    GRanges(seqnames = x$chr,
                  ranges = IRanges(
                      start = as.numeric(x$start),
                      end = as.numeric(x$end)))})

## want regions that overlap input and ChIP-peaks
input_peaks = findOverlaps(peaksgr[["Input"]],peaksgr[["Peaks"]])
input_subset = peaksgr[["Input"]][queryHits(input_peaks)]
peaks_subset = peaksgr[["Peaks"]][subjectHits(input_peaks)]

readfiles = readfiles[grep("K562",readfiles)]
readfiles = readfiles[grep("sort",readfiles)]
readfiles = readfiles[grep("bai",readfiles,invert = TRUE)]


reads = GenomicAlignments::readGAlignments(readfiles,param = NULL)
reads = as(reads,"GRanges")

exo = ExoData(reads = reads,mc.cores = 20,ntimes = 1e3)

library(ggplot2)
library(hexbin)
library(scales)

figs = "figs/NAR_review"

ov = findOverlaps(exo, input_subset)


source("~/Desktop/Docs/Code/ChIPexoQual/R/base_summaryStats.R")

ntimes = 1000
nregions = 1000

exoDT = exo %>% as.data.frame %>% as.tbl %>%
    select(uniquePos,width,depth) %>%
    as.data.table


param_runs = list()

param_runs[["all"]] = do.call(rbind,
                              mclapply(seq_len(ntimes),
                                       calculateParamDist,
                                       exoDT,nregions,mc.cores = 20))
idx = ov %>% queryHits %>% unique                                       


param_runs[["no_artifact"]] =
    do.call(rbind,
            mclapply(seq_len(ntimes),
                     calculateParamDist,
                     exoDT[-idx],nregions,mc.cores = 20))

param_runs[["artifact"]] =
    do.call(rbind,
            mclapply(seq_len(ntimes),
                     calculateParamDist,
                     exoDT[idx],nregions,mc.cores = 20))


ll = mapply(function(x,y){
    x = x[,set := y]
    x},param_runs,names(param_runs),SIMPLIFY = FALSE)
ll = rbindlist(ll)

ll[,set := factor(set,levels = c("all","no_artifact","artifact"))] 

library(ggplot2)

theme_set(theme_bw())

pdf(file.path(figs,"GR_artifact.pdf"))
ggplot(ll[term == "uniquePos"],
       aes(set,estimate,colour = set))+geom_boxplot()+
    ylab(expression(beta[1]))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")+ylim(0,100)
ggplot(ll[term == "uniquePos"],
       aes(set,estimate,colour = set))+geom_boxplot()+
    ylab(expression(beta[1]))+scale_y_log10()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")
ggplot(ll[term == "width"],
       aes(set,-estimate,colour = set))+geom_boxplot()+
    ylab(expression(beta[2]))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "top")+scale_color_brewer(palette = "Set1",name = "")+ylim(-.1,5)
dev.off()

