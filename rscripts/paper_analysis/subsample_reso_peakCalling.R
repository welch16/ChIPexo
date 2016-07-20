
rm(list = ls())

library(mosaics)
library(data.table)
library(ggplot2)
library(parallel)
library(GenomicRanges)
library(data.table)

## param
FDR = 0.05
bin_size = 150
thres = 10

indir = "/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity"

files = list.files(indir,recursive = TRUE,full.names = TRUE)
files = files[grep("bins",files)]

input = files[grep("INPUT",files)]
files = files[grep("INPUT",files,invert = TRUE)]

exo = files[grep("ChIPseq",files,invert = TRUE)]
pet = files[grep("PET",files)]
set = files[grep("SET",files)]

get_peak = function(fit){
  mosaicsPeak(fit,signalModel = ifelse(fit@bic2S <= fit@bic1S,
                   "2S","1S"),FDR = FDR,maxgap = bin_size,
              thres = thres)}

## ChIPexo
mappability = "/p/keles/genome_data/EColi_U00096.2/mappability/bin/fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt"
gc_content = "/p/keles/genome_data/EColi_U00096.2/GC/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt"
naked_dna = "/p/keles/genome_data/EColi_U00096.2/N/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"

exobins = mclapply(exo,function(x)
      readBins(type = c("chip","M","GC","N"),
         fileName = c(x,
           mappability,
           gc_content,
           naked_dna)),mc.cores = 4)

exofit = mclapply(exobins,mosaicsFit,mc.cores = 4)

par(mfrow = c(2,2))
lapply(exofit,plot)

## ChIPseq PET
petbins = mclapply(pet,function(x)
  readBins(type = c("chip","input"),
           fileName = c(x,input)),mc.cores = 2)

petfit = mclapply(petbins,mosaicsFit,mc.cores = 2)

## ChIPseq SET
setbins = mclapply(set,function(x)
  readBins(type = c("chip","input"),
           fileName = c(x,input)),mc.cores = 2)
setfit = mclapply(setbins,mosaicsFit,mc.cores = 2)

par(mfrow = c(2,2))
lapply(petfit,plot)
lapply(setfit,plot)

exopeak = mclapply(exofit,get_peak,mc.cores = 4)
names(exopeak) = basename(exo)

petpeak = mclapply(petfit,get_peak,mc.cores = 2)
names(petpeak) = basename(pet)

setpeak = mclapply(setfit,get_peak,mc.cores = 2)
names(setpeak) = basename(set)

peaks = list(exo = exopeak,pet = petpeak,set = setpeak)

save(peaks,file = file.path(indir,"subsample_reso_peaks.RData"))
