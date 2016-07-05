
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(data.table)
library(parallel)

bin_size <- 100
frag_len <- 150

dr <- "/p/keles/ChIPexo/volume6/imbalance"

files <- list.files(dr,full.names = TRUE,recursive = TRUE)
files <- files[grep("chr",files,invert = TRUE)]

addfiles <- files[grep("addf",files)]
files <- files[grep("addf",files,invert = TRUE)]
names(addfiles) <- c("GC","M","N")

input <- files[grep("input",tolower(files))]
files <- files[grep("input",tolower(files),invert = TRUE)]


fit_mosaics_OS <- function(binfile,extra,bgEst)
{
  type <- c("chip","M","GC","N")
  bins <- readBins(type = type, fileName = c(binfile,extra[type[-1]]),
                   parallel = TRUE,nCore = 24) 
  fit <- mosaicsFit(bins,bgEst=bgEst,parallel = TRUE,nCore = 24)
  return(fit)
}

fit_mosaics_IO <- function(binfile,extra,bgEst = "rMOM")
{
  type <- c("chip","input")
  bins <- readBins(type = type, fileName = c(binfile,extra),
                   parallel = TRUE,nCore = 24)
  fit <- mosaicsFit(bins,bgEst = bgEst,parallel = TRUE,nCore = 24 )
  return(fit)
}


## CTCF - exo
ctcf_files <- files[grep("ctcf",tolower(files))]

fit_CTCF <- list()
fit_CTCF[["exo"]] <- fit_mosaics_OS(ctcf_files[1],
   addfiles,bgEst = "rMOM") ## even work with rMOM which is the suggestion
## when for well sequenced experiments
fit_CTCF[["seq1"]] <- fit_mosaics_IO(ctcf_files[2],input[1],bgEst = "matchLow")
fit_CTCF[["seq2"]] <- fit_mosaics_IO(ctcf_files[2],input[1],bgEst = "rMOM")


pdf(file = "figs/ctcf_gof.pdf")
u <- lapply(fit_CTCF,plot)
dev.off()


call_peaks <- function(fit,thres,fdr)
{
  mod <- ifelse(fit@bic1S < fit@bic2S ,"1S","2S")
  mp <- mosaicsPeak(fit,signalModel = mod,thres = 100,FDR = fdr)
  peaks <- data.table(print(mp))
  return(peaks)
}

thres <- 100
fdr <- 5

peaks <- lapply(fit_CTCF,call_peaks,thres , fdr)

save(peaks,file = "data/CTCF_peaks.RData")

#fit_ER <- lapply(files[grep("ERR",files)],fit_mosaics,addfiles,"rMOM")

## fit_TBP <- lapply(files[grep("TBP_K562",files)],fit_mosaics,addfiles,"matchLow")

## #peaks_ER <- lapply(fit_ER,call_peaks,thres,fdr)
## peaks_TBP <- lapply(fit_TBP,call_peaks,thres,fdr)

## #save(peaks_ER,file = "data/ChIPexo_ER_peaks.RData")
## save(peaks_TBP,file = "data/ChIPexo_TBP_peaks.RData")
