rm(list = ls())

library(ChIPexoQual)
library(magrittr)
library(readr)
library(dplyr)
library(data.table)

dr = "/p/keles/ChIPexo/volume4"

files = list.files("data/figures/fig5",full.names  = TRUE,pattern = "TBP")
files = files[grep("0M.",files,invert = TRUE)]

scorefiles = files[grep("scores",files)]
statfiles = files[grep("stat",files)] 

names(scorefiles) = c(paste0("Nexus",seq_len(2)),
                      paste0("Exo",seq_len(3)))

names(statfiles) = names(scorefiles)

stats = lapply(statfiles,read_tsv)
scores = lapply(scorefiles,read_tsv)

peakfiles = list.files(dr,full.names = TRUE,recursive = TRUE,pattern = "peaks")
peakfiles = peakfiles[grep("venters",peakfiles)]
peakfiles = peakfiles[grep("encode",peakfiles,invert = TRUE)]

peaks = mclapply(peakfiles,read_delim,delim = " ",col_names = FALSE)
peaks = lapply(peaks,function(x){
    x = x %>% select(X1,X2,X3)
    setnames(x,names(x),c("seqnames","start","end"))
    x = as.data.table(x)
  return(ChIPUtils::dt2gr(x))})

## peak columns:
## chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC 

## join all peak regions together
all_peaks = Reduce(c,peaks) %>% reduce

exo = stats %>% lapply(function(x){
    gr = GRanges(seqnames = x$seqnames,
            ranges = IRanges(
                start = x$start,
                end = x$end))
    mcols(gr) = x[,-c(1:5)] %>% DataFrame
    gr})

exo_peaks = mclapply(exo,subsetByOverlaps,all_peaks)

readlength = exo %>% sapply(function(x)x %>% width %>% min)

## common sense filter

## remove chrM
exo_peaks = exo_peaks %>% mclapply(function(x)x[as.character(seqnames(x)) != "chrM"])

## width analysis

exo_peaks %>% lapply(function(x)width(x) %>% summary)


rl = median(readlength)

K = 3
exo_peaks = exo_peaks %>% lapply(function(x)x[width(x) >= K * rl])

## load fimo files
fimofiles = list.files(dr,pattern = "fimo",recursive = TRUE,full.names = TRUE)
fimofiles = fimofiles[grep("txt",fimofiles)]
fimofiles = fimofiles[grep("fimo_peaks",fimofiles)]
fimofiles = fimofiles[grep("FOXA1",fimofiles,invert = TRUE)]
fimo = lapply(fimofiles,read_delim,delim = "\t")
names(fimo) = c(paste0("Nexus",1:2),paste0("Exo",1:3))

fimo = fimo[names(exo)]

library(ggplot2)
library(tidyr)

## scores boxplots
fimo = mapply(function(x,y)x %>%
                           mutate(replicate = y),fimo,names(fimo),SIMPLIFY = FALSE)


topK = c(50,100,250,500,1000,2000,4000,8000)


dt_list = fimo %>% lapply(function(x){
    lapply(topK,function(z)top_n(x,z,score) %>%
                           mutate(K = z)) %>% bind_rows}) %>% bind_rows

write_tsv(dt_list,"data/figures/fig4/fig4_FIMO_scores_TBP.tsv")





strand_imbalance_fimo <- function(exo,fimo,repl)
{

    exo_tbl = exo %>% as.data.frame %>% as.tbl %>%
        mutate(match = paste0(seqnames,":",start,"-",end)) %>%
        select(match , everything())
    fimo = fimo %>% rename(match = `sequence name`)

    fimo_summary = fimo %>% group_by(match) %>%
        summarize(
            nmotif = n(),
            minPval = min(`p-value`),
            maxScore = max(score)
        )

    exo_tbl %>% left_join(fimo_summary,by ="match") %>%
        mutate(nmotif = ifelse(is.na(nmotif),0,nmotif),repl)

  
}


imbalance = mapply(strand_imbalance_fimo,exo_peaks,fimo,names(fimo),SIMPLIFY = FALSE) %>% bind_rows


write_tsv(imbalance, "data/figures/fig4/fig4_imbalance_hist_TBP.tsv")
