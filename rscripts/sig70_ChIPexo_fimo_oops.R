
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(RColorBrewer)

fimo_dr <- "/p/keles/ChIPexo/volume6/K12/meme_sig70"
files <- list.files(fimo_dr,pattern = "fimo.txt",
                    full.names = TRUE,recursive = TRUE)

fimo <- lapply(files,fread)
fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("pattern","peakID",
                        "start","end","strand","score","pval",
                        "qval","sequence"))
  return(x)})
names(fimo) <- sapply(strsplit(basename(dirname(files)),"_"),
       function(x)x[1])

