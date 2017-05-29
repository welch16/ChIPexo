
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


## ER - exo
tbp_files <- files[grep("TBP",files)]

fit_TBP <- list()
fit_TBP[["exo1"]] <- fit_mosaics_OS(tbp_files[3],
   addfiles,bgEst = "matchLow")
fit_TBP[["exo2"]] <- fit_mosaics_OS(tbp_files[4],
   addfiles,bgEst = "matchLow")
fit_TBP[["exo3"]] <- fit_mosaics_OS(tbp_files[5],
   addfiles,bgEst = "matchLow")
fit_TBP[["seq1"]] <- fit_mosaics_IO(tbp_files[1],
  input[2],bgEst = "matchLow")                                    
fit_TBP[["seq2"]] <- fit_mosaics_IO(tbp_files[2],
  input[3],bgEst = "matchLow")                                    

pdf(file = "figs/TBP_gof.pdf")
u <- lapply(fit_TBP,plot)
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

peaks <- lapply(fit_TBP,call_peaks,thres , fdr)

save(peaks,file = "data/TBP_peaks.RData")

