
rm(list = ls())

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(parallel)
library(GenomicAlignments)
library(readr)

devtools:::load_all("~/Desktop/Docs/Code/mosaics")

figs_dir = "/p/keles/ChIPexo/volume3/ChIPexo/figs/NAR_review/TBPexo_peaks"

dir.create(figs_dir,showWarnings = FALSE)

work_dir = "/p/keles/ChIPexo/volume4/venters_data/sortbam"

infiles = list.files(work_dir,full.names = TRUE,recursive = TRUE)
infiles = infiles[grep("TBP",infiles)]
infiles = infiles[grep("bai",infiles,invert = TRUE)]

reads = mclapply(infiles,ChIPUtils:::create_reads,mc.cores =8)

chrom.sizes = read_delim("/p/keles/SOFTWARE/hg19.chrom.sizes",delim = "\t",col_names =FALSE)
chrom.sizes = GRanges(seqnames = chrom.sizes$X1,
                      ranges = IRanges(start = 1 , end = chrom.sizes$X2))

bins = mclapply(reads,function(x){
    ChIPUtils:::create_bins(200,reads = x,chrom = chrom.sizes, frag_len = 200)},mc.cores = 8)

binsDF = lapply(bins,function(x)tibble(chr = as.character(seqnames(x)),start = as.integer(start(x)),
                                       counts = as.integer(x$tagCounts)))

depth = lapply(reads,function(x)x@nReads)

write_bins <- function(bin,depth,file)
{
    line = paste(sep = " ","#","sequencing depth:",depth)
    write(line,file = file,append = TRUE)

    write_delim(bin,path = file,delim = "\t",append = TRUE)

}

u = mcmapply(write_bins,binsDF,depth,file.path(work_dir,"bins",
         paste0(basename(infiles),"_fragL",200,"_bin",200,".txt")),
         mc.cores = 8)

## first step, build data

read_bins <- function(file)
{
    Mfile  = file.path(work_dir,"extra","all_map_fragL200_bin200.txt")
    GCfile = file.path(work_dir,"extra","all_GC_fragL200_bin200.txt")
    Nfile = file.path(work_dir,"extra","all_N_fragL200_bin200.txt")
    readBins(c("chip","M","GC","N"),
             c(file,Mfile,GCfile,Nfile))
}

binfiles = list.files(file.path(work_dir,"bins"),full.names = TRUE)

bins = mclapply(binfiles,read_bins,mc.cores = 8)

