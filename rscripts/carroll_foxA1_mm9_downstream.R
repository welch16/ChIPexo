
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(dpeak)
library(data.table)

indir <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(indir)

files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]

outdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"

### construct bins
fragLen <- 200
binSize <- 200

A <- mclapply(file.path(indir,files),
  constructBins,fileFormat = "bam",
              outfileLoc = file.path(outdir,"bins"),
              fragLen = fragLen ,binSize = binSize,mc.cores = 3)

join_bins <- function(dr)
{
  mc <- detectCores()
  files <- list.files(dr)
  sqnms <- sapply(strsplit(files,"_",fixed = TRUE),function(x)x[1])
  bins <- mclapply(file.path(dr,files),read.table,mc.cores = mc)
  bins <- lapply(bins,data.table)
  bins <- mcmapply(function(sqn,bin){
    setnames(bin,names(bin),c("start","val"))
    bin[,seqnames := sqn]
    setcolorder(bin,c("seqnames","start","val"))
    return(bin)
  },sqnms,bins,SIMPLIFY = FALSE,mc.cores = mc)
  bins <- do.call(rbind,bins)
  tt <- tempfile(pattern = "bin",fileext = ".txt")
  write.table(bins,file = tt ,col.names = FALSE,quote = FALSE,row.names = FALSE,sep = "\t")
  return(tt)
}

A <- lapply(file.path(outdir,c("GC","N","mappability")),join_bins)
names(A) <- c("GC","N","map")

bindir <- file.path(outdir,"bins")
binfiles <- list.files(bindir)

bins <- lapply(file.path(bindir,binfiles),function(x)readBins(type = c("chip","M","GC","N"),
  fileName = c(x,A[["map"]],A[["GC"]],A[["N"]]),
  parallel = TRUE,nCore = detectCores()))

                                  
fits <- lapply(bins,mosaicsFit,parallel = TRUE,nCore = detectCores())
save(fits,file = "./fits.RData")

pdf(file = file.path(outdir,"mosaics_gof","GOF_plots.pdf"))
lapply(fits,plot)
dev.off()

load("/fits.RData")

call_peaks <- function(fit,fdr,thres,mg)
{
  mp <- mosaicsPeak(fit,signalModel = "2S",thres = thres,maxgap = mg)
  peaks <- data.table(print(mp))
  
  return(peaks)
}

fdr <- 0.05
thres <- 100
mg <- 200

peaks <- lapply(fits,call_peaks,fdr,thres,mg)

## peak distribution by replicate and chr
## > sapply(peaks,function(x)x[,summary(chrID)])
##       [,1] [,2] [,3]
## chr1   490  872  267
## chr10  360  743  189
## chr11  453  779  221
## chr12  288  486  158
## chr13  336  663  178
## chr14  286  511  134
## chr15  319  563  174
## chr16  264  499  133
## chr17  239  440  127
## chr18  251  487  127
## chr19  259  515  112
## chr2   478  872  242
## chr3   403  755  229
## chr4   447  821  223
## chr5   473  816  232
## chr6   419  818  211
## chr7   318  581  140
## chr8   382  668  216
## chr9   375  674  196
## chrM     2    3    2
## chrX    72  143   33
## chrY    34   40   40
###---------------------
##       6948 12749  3584


mapply(write.table,peaks,file = file.path(outdir,"peaks",gsub(".sort.bam", "_peaks.txt",files)),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

## chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC  





