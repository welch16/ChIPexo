
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(devtools)
library(dplyr)

load("data/sig70_ChIPSeq_PE_unified_peaks.RData") ## PE ChIP-Seq peaks

dr <- "/p/keles/ChIPexo/volume6/K12/meme_sig70"

fimo_files <- list.files(dr,recursive = TRUE,pattern = "fimo.txt",
                   full.names = TRUE)
fimo <- lapply(fimo_files,fread)


fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("pattern","seqID",
                        "start","end","strand",
                        "score","pvalue","qvalue","sequence"))
  return(x)})

names(fimo) <- basename(dirname(fimo_files))

fimo <- lapply(fimo,function(x){
  seq <- t(x[,strsplit(gsub("U00096:","",seqID),"_",fixed = TRUE)])
  x[,regionStart := as.numeric(seq[,1])]
  x[,regionEnd := as.numeric(seq[,2])]
  return(x)
})


edsn <- c("931","933","935","937","1311","1314","1317","1320")

fimo <- lapply(edsn,function(x,fimo){
  out <- fimo[grep(x,names(fimo))]
  pat <- lapply(out,function(z)z[,unique(pattern)])
  out[[2]][,pattern := plyr::mapvalues(pattern,from = pat[[2]],
              to = pat[[2]] + max(pat[[1]]))]
  return(do.call(rbind,out))  
},fimo)
names(fimo) <- edsn

## get_peaks <- function(x){
##   seq <- x[,seqID]
##   seq <- strsplit(seq,":",fixed = TRUE)
##   sqnms <- sapply(seq,function(z)z[1])
##   seq <- sapply(seq,function(z)z[2])
##   seq <- strsplit(seq,"_",fixed = TRUE)
##   st <- as.numeric(sapply(seq,function(z)z[1]))
##   en <- as.numeric(sapply(seq,function(z)z[2]))
##   out <- sort(reduce(GRanges(seqnames = sqnms,
##                         ranges = IRanges(start = st,
##                           end = en))))
##   return(out)                  
## }


## peaks <- lapply(fimo,get_peaks)
