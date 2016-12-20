
rm(list = ls())

library(ChIPexoQual)
library(dplyr)
library(magrittr)
library(readr)
library(parallel)
library(GenomicRanges)
library(tidyr)
library(ggplot2)

devtools:::load_all("~/Desktop/Docs/Segvis")

## read fimo output
fimo_dir = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(fimo_dir,full.names=  TRUE,recursive = TRUE)
files = files[grep("FOXA1",files)]
files = files[grep("txt",files)]
files = files[grep("fasta",files)]
peaksMotif = lapply(files,read_delim,delim = "\t")

## read aligned reads
read_dir = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
rfiles = list.files(read_dir,full.names = TRUE,recursive = TRUE)
rfiles = rfiles[grep("sort",rfiles)]
rfiles = rfiles[grep("txt",rfiles,invert = TRUE)]
rfiles = rfiles[grep("bai",rfiles,invert = TRUE)]

## generate GRanges for peaks
load("data/FoxA1_ChIPexo_peaks.RData")  # peaks
peaks = as.tbl(peaks)

peak_gr <- function(repl , peaks)
{
    out = peaks %>% filter(rep == repl)
    GRanges(seqnames = out$chrID,
            ranges = IRanges(
                start = out$peakStart,
                end = out$peakStop))
}

regions  = lapply(c("rep3","rep1","rep2"),peak_gr,peaks)
                    
motifs = lapply(peaksMotif,function(x){
    seq = x %>% select(contains("name")) %>% select(contains("sequen"))
    names(seq) = "seq"
    seq = seq %>% separate(seq,c("seqnames","range"),sep = ":") %>%
        separate(range,c("start","end"),sep = "-")
    gr = GRanges(seqnames = seq$seqnames,
            ranges = IRanges(
                start = as.numeric(seq$start),
                end = as.numeric(seq$end)))
    resize(gr,1,fix = "center")
})

which.motif <- function(peak,motif)
{
    ov = findOverlaps(peak,motif)
   peak[queryHits(ov)]
}

regions.motif = mapply(which.motif,regions,motifs,SIMPLIFY = FALSE)

reads = mclapply(rfiles,readGAlignments,param = NULL,mc.cores = 3)
reads = mclapply(reads,as,"GRanges",mc.cores = 3)

calculate_strand_imbalance <- function(region.m ,read)
{
    fwd = subset(read,as.character(strand(read)) == "+")
    f = countOverlaps(region.m,fwd)
    all = countOverlaps(region.m , read)
    f / all
}

## calculate stats
fsr = mcmapply(calculate_strand_imbalance,regions.motif,reads,SIMPLIFY = FALSE,mc.cores = 3)
depths = mcmapply(countOverlaps,regions.motif,reads,SIMPLIFY = FALSE,mc.cores = 3)

## plots with segvis
segvislist = mapply(SegvizData,regions.motif,rfiles,
                    MoreArgs = list(frag_len = 36),SIMPLIFY = FALSE)


mytitle <- function(gr)
{
    paste0(as.character(seqnames(gr)),":", prettyNum(start(gr),big.mark = ","),"-",
           prettyNum(end(gr),big.mark = ","))
}

theme_set(theme_bw())

plot_regions <- function(reg,segvis,motifs,fsr,dept,name)
{
    mm = findOverlaps(segvis[reg],motifs)
    mm = motifs[subjectHits(mm)]
    plot_region(segvis,reg,type = "both",nameFiles = name)+
        theme(legend.position = "top")+
        labs(title = mytitle(segvis[reg]),
             subtitle = paste("FSR:",round(fsr[reg],4),"and Depth:",prettyNum(dept[reg],big.mark = ",")))+
        geom_vline(xintercept = mid(ranges(mm)),linetype = 2)       
}

## demo = mclapply(seq_len(120),plot_regions,segvislist[[1]],
##                 motifs[[1]],fsr[[1]],depths[[1]],"Rep-3",mc.cores = 12)

rep3 = mclapply(seq_len(length(segvislist[[1]])),plot_regions,segvislist[[1]],
                motifs[[1]],fsr[[1]],depths[[1]],"Rep-3",mc.cores = 12)

rep1 = mclapply(seq_len(length(segvislist[[2]])),plot_regions,segvislist[[2]],
                motifs[[2]],fsr[[2]],depths[[2]],"Rep-1",mc.cores = 12)

rep2 = mclapply(seq_len(length(segvislist[[3]])),plot_regions,segvislist[[3]],
                motifs[[3]],fsr[[3]],depths[[3]],"Rep-2",mc.cores = 12)


figs_dir = "figs/NAR_review/FSR_peaks"

## pdf(file.path(figs_dir,"FoxA1_rep3_demo.pdf"),width = 7,height = 5)
## u = lapply(demo,print)
## dev.off()


pdf(file.path(figs_dir,"FoxA1_rep3.pdf"),width = 7,height = 5)
u = lapply(rep3,print)
dev.off()

pdf(file.path(figs_dir,"FoxA1_rep1.pdf"),width = 7,height = 5)
u = lapply(rep1,print)
dev.off()

pdf(file.path(figs_dir,"FoxA1_rep2.pdf"),width = 7,height = 5)
u = lapply(rep2,print)
dev.off()




