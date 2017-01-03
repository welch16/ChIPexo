
rm(list = ls())

library(magrittr)
library(dplyr)
library(readr)
library(data.table)

dr = "/p/keles/ChIPexo/volume4/meijsing_data"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("xls",files)]

load("data/ChIPexo_QC_runs/meijsing_GR_K562_Rep1.RData")

peaks = lapply(files,read_delim,delim = "\t",comment = "# ")
peaks[[1]] = peaks[[1]][-1,]
peaks[[3]] = peaks[[3]][-1,]

names(peaks) = c("Input","Negative","Peaks")

library(GenomicRanges)

stats = ext_stats$stats

statsgr = stats[,GRanges(seqnames = seqnames,
                       ranges = IRanges(
                           start = start,
                           end = end))]

peaksgr = lapply(peaks,function(x){
    GRanges(seqnames = x$chr,
                  ranges = IRanges(
                      start = as.numeric(x$start),
                      end = as.numeric(x$end)))})

## want regions that overlap input and ChIP-peaks
input_peaks = findOverlaps(peaksgr[["Input"]],peaksgr[["Peaks"]])
input_subset = peaksgr[["Input"]][queryHits(input_peaks)]
peaks_subset = peaksgr[["Peaks"]][subjectHits(input_peaks)]

stats_input = findOverlaps(statsgr,input_subset)
stats_peaks = findOverlaps(statsgr,peaks_subset)

stats1 = stats[queryHits(stats_input)]
stats2 = stats[queryHits(stats_peaks)]

stats1[,which.max(depth)]

library(GenomicAlignments)

devtools:::load_all("~/Desktop/Docs/Segvis/")

gr = stats1[,GRanges(seqnames = seqnames,ranges = IRanges(start = start,end = end))]

## nrow(stats1) / nrow(stats) * 100
## [1] 0.6351425

library(ggplot2)
library(hexbin)
library(scales)

figs = "figs/NAR_review"


## pdf(file = file.path(figs,"ARC_vURC.pdf"))
## stats1 %>% ggplot(aes(ave_reads,cover_rate))+stat_binhex(bins = 100)+
##     scale_fill_gradientn(colours = r, trans = "log10")+xlim(0,10)
## dev.off()


ss = SegvizData(regions = gr[33748],files = "/p/keles/ChIPexo/volume4/meijsing_data/K562_GR_chip-exo.sort.bam",frag_len = 36,mc.cores = 12)


pdf(file.path(figs,"GR_regions.pdf"),height = 5)
print(plot_region(ss,1,nameFiles = "GR_K562",type = "both",normalize =  FALSE))
dev.off()
