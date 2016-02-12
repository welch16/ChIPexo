
rm(list = ls())

library(dpeak)
library(GenomicAlignments)
library(data.table)
library(parallel)

peakdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/peaks"
peakfiles <- list.files(peakdir)

readdir <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(readdir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]

nCore <- detectCores()

maxG <- 5

dpeak_sites <- function(peaksfile,readsfile,outfile)
{
  dpeak <- dpeakRead(peakfile = peaksfile,readfile = readsfile,fileFormat ="bam",PET = FALSE,nCore = nCore)
  fit <- dpeakFit(dpeak, maxComp = maxG,nCore = nCore)
  tt <- tempfile(fileext = "bed")
  export(fit,type = ".bed",filename = tt)
  dt <- data.table(read.table(tt,skip = 1))
  setnames(dt,names(dt),c("chrID","start","end","peakID","strength"))
  write.table(dt,file = outfile,quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
}

outfiles <- gsub("_peaks","_sites",peakfiles)
outdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/sites"

mapply(dpeak_sites,file.path(peakdir,peakfiles),
       file.path(readdir,files),
       file.path(outdir,outfiles))
