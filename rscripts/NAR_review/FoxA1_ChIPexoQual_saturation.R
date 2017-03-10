
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

dr = "/p/keles/ChIPexo/volume4"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("carro",files)]
files = files[grep("sort",files)]
files = files[grep("mouse",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("txt",files,invert = TRUE)]

library(parallel)
set.seed(12345)

## reads = readGAlignments(files[5],param = NULL)
## reads = as(reads,"GRanges")
## nreads = reads %>% length

inc = 4e6
maxdepth = 2e7

sample_reads <- function(reads,inc,maxdepth)
{

    size = floor(length(reads) / inc)
    idd = sample( size,length(reads),replace = TRUE)

    v = seq_len(size) * inc
    maxgroups = which.max(v[v <= maxdepth])

    out  = lapply(seq_len(maxgroups),function(i){
        idx = which(idd <= i)
        reads[idx]})
    names(out) = v[seq_along(out)] %>% as.integer %>% prettyNum(big.mark = ",")
    out
    
}

## sr = sample_reads(reads,inc,maxdepth)

reads = files %>% lapply(readGAlignments,param = NULL)
reads = reads %>% lapply(as,"GRanges")
names(reads) = paste("Rep",c(3,1,2))


samp_reads = lapply(reads,sample_reads,inc,maxdepth)

samp_reads = samp_reads %>% unlist

base = lapply(reads,function(x)ExoData(reads = x,mc.cores = 24))

exo = lapply(samp_reads,function(x)ExoData(reads = x, mc.cores = 24))


get_overlap <- function(dataset,base,exo)
{
    my_base = strsplit(dataset,"\\.")[[1]][1]
    ov = findOverlaps(
        exo[[dataset]],
        base[[my_base]]
    )
    ov 
}

overlaps = lapply(names(samp_reads),get_overlap,base,exo)
names(overlaps) = names(samp_reads)


overlap_proportion <- function(ov)
{
    a = subjectHits(ov) %>% unique %>% length
    a / subjectLength(ov)
}

tbl = tibble(nm = names(overlaps),
             prop = overlaps %>% sapply(overlap_proportion)
             )

tbl = tbl %>% separate(nm,into = c("Dataset","samp"),sep = "\\.") %>%
    mutate(sample = gsub(",","",samp) %>% as.numeric)


pdf("figs/NAR_review/threshold/FoxA1_saturation.pdf")
u = tbl %>% ggplot(aes(sample,prop,colour = Dataset))+geom_line(size = 1.5)+
    scale_x_continuous(breaks = tbl$sample %>% unique ,
                       labels = tbl$sample %>% unique %>% as.integer %>% prettyNum(big.mark = ","))+
    theme_bw()+theme(legend.position = "top",axis.text.x = element_text(angle = 15,hjust =  1))+
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0,.25,.5,.75,.9,.95,1),
                       limits = c(0,1))+
    scale_color_brewer(palette = "Set1")+xlab("Number of aligned reads used")+
    ylab("Percentage of regions shared with whole dataset")
u = print(u)
dev.off()                       


pdf("figs/NAR_review/threshold/FoxA1_Depth_UniquePosRatio_histogram.pdf")



dev.off()
