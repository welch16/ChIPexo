
rm(list = ls())

library(parallel)
library(data.table)
library(GenomicAlignments)

## library(devtools)
## load_all("~/Desktop/Docs/Code/dpeak")
library(dpeak)

mc <- 12

peak_dir <- "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_SET/peaks"
peak_files <- list.files(peak_dir)


read_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET/"
read_files <- list.files(read_dir)
read_files <- read_files[grep("filter.bam",read_files)]
read_files <- read_files[grep("bai",read_files,invert = TRUE)]
read_files <- read_files[-1]


fragLen <- 150

dpeak_read_wrap <- function(peak,read_files,peak_dir,read_dir,fragLen)
{
  id <- strsplit(peak,"_")[[1]][1]
  message(peak)
  readfile <- read_files[grep(id,read_files)]
  out <- dpeakRead(peakfile = file.path(peak_dir,peak),
                   readfile = file.path(read_dir,readfile),
                   fileFormat = "bam",fragLen = fragLen)
  return(out)
}


dpeaks <- lapply(peak_files[1],dpeak_read_wrap,read_files,peak_dir,read_dir,fragLen)

maxComp <- 5
fits <- lapply(dpeaks,dpeakFit,maxComp = maxComp,nCore = mc)

out_dir <- "/p/keles/ChIPexo/volume6/results/dpeak/Landick/ChIPseq_SET"

export_wrap <- function(peak_file,fit,out_dir)
{
  id <- strsplit(peak_file,"_")[[1]][1]
  export(fit,type = "bed",filename = file.path(out_dir,paste0(id,"_binding_sites_common_maxComp",maxComp ,".bed")))

}

out <- mapply(export_wrap,peak_files,fits,MoreArgs = list(out_dir))
