rm(list = ls())

library(mosaics)
library(data.table)
library(GenomicAlignments)
library(dpeak)

######################################################################################

## Condition table

source("R/base_edsn.R")

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab_old)
names(char) <- what

#############1#########################################################################

## Initial parameters

tf <- "Sig70"
cond <- "aerobic"
bs <- 150
fl <- 150
fdr <- .1
thresh <- 10
mc <- 24
maxComp <- 1

flag_bins <- FALSE
flag_peaks <- FALSE
flag_binding <- TRUE
flag_input <- FALSE

char <- lapply(char,function(x)x[ip == tf & condition == cond])

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

base_dir <- "/p/keles/ChIPexo/volume6/condition"
folder <- paste(tf,cond,sep = "_")

base_dir <- file.path(base_dir,folder)
check_create(base_dir)

##################################h####################################################

## Input data bins for ChIPseq (PET and SET)

inputs_dir <- "/p/keles/ChIPexo/volume6/inputs"
inputs_dirs <- lapply(what[-1],function(x)file.path(inputs_dir,x))

chip_dirs <- list()
chip_dirs[["exo"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
chip_dirs[["pet"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
chip_dirs[["set"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"


if(flag_input){
  lapply(inputs_dirs,check_create)
  constructBins(file.path(chip_dirs[["pet"]],"edsn1416_042814_qc.filter.bam"),
    fileFormat = "bam",outfileLoc = inputs_dirs[[1]],
    PET = TRUE,fragLen = fl,binSize = bs)
  constructBins(file.path(chip_dirs[["set"]],"edsn101_042814_qc.sorted.bam"),
    fileFormat = "bam",outfileLoc = inputs_dirs[[2]],
    PET = FALSE,fragLen = fl,binSize = bs)                              
}

######################################################################################

### build_bins

bin_dir <- file.path(base_dir,"bins")
check_create(bin_dir)

bin_dirs <- lapply(what,function(x)file.path(bin_dir,x))
lapply(bin_dirs,check_create)

sample_read_files <- function(in_dir,dt,pet)
{
  files <- list.files(in_dir)
  edsn <- dt[,(edsn)]
  files <- files[sapply(edsn,function(x)grep(x,files))]
  if(pet){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("bai",files,invert = TRUE)]
  files <- files[grep("run",files,invert = TRUE)]
  return(files)
}

create_bins <- function(in_dir,out_dir,dt,pet,bs,fl)
{
  files <- sample_read_files(in_dir,dt,pet)
  lapply(file.path(in_dir,files),constructBins,
         fileFormat = "bam",
         outfileLoc = out_dir,
         PET = pet,
         fragLen = fl,
         binSize = bs)                  
}

if(flag_bins){
  mapply(create_bins,chip_dirs,bin_dirs,char,list(FALSE,TRUE,FALSE),
    MoreArgs = list(bs,fl),SIMPLIFY = FALSE)
}

######################################################################################

## call peaks

peak_dir <- file.path(base_dir,"peaks")
check_create(peak_dir)

peak_dirs <- lapply(what,function(x)file.path(peak_dir,x))
lapply(peak_dirs,check_create)

peak_dirs <- file.path(peak_dirs,paste0("FDR",fdr*100))
lapply(peak_dirs,check_create)

extra <- list()
extra[["exo"]] <- file.path("/p/keles/genome_data/EColi_U00096.2",
  c("mappability","GC","N"),"bin",
  c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))
extra[["pet"]] <- "/p/keles/ChIPexo/volume6/inputs/pet/edsn1416_042814_qc.filter.bam_bin150.txt"
extra[["set"]] <- "/p/keles/ChIPexo/volume6/inputs/set/edsn101_042814_qc.sorted.bam_fragL150_bin150.txt"

call_peaks <- function(file,opt,extra,fdr,bs,thresh)
{
##  message(file)
  bins <- readBins(opt, c(file,extra))
  fit <- mosaicsFit(bins)
  if(fit@bic2S <= fit@bic1S){
    peak <- mosaicsPeak( fit, signalModel="2S", FDR= fdr, maxgap=bs, thres=thresh )
  }else{
    peak <- mosaicsPeak( fit, signalModel="1S", FDR= fdr, maxgap=bs, thres=thresh )
  }
  peak <- data.table(peak@peakList)
  return(peak)
}

mosaics_peaks_wrap <- function(in_dir,out_dir,what,extra,fdr,thresh,bs,mc)
{
  if(what == "exo"){
    opt <- c("chip","M","GC","N")
  }else{
    opt <- c("chip","input")
  }

  files <- file.path(in_dir,list.files(in_dir))
  peaks <- mclapply(files,call_peaks,opt,extra,fdr,bs,thresh,mc.cores = mc)

  files <- list.files(in_dir)
  files <- sapply(strsplit(files,".",fixed = TRUE),function(x)x[1])
  files <- file.path(out_dir,paste0(files,"_peaks.txt"))

  out <- mcmapply(write.table,peaks,files,MoreArgs = list(
    sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE),SIMPLIFY = FALSE,mc.cores = mc)  

  return(out)

}

if(flag_peaks){
  z <- mapply(mosaics_peaks_wrap,bin_dirs,peak_dirs,what,extra,
    MoreArgs = list(fdr,thresh,bs,mc),SIMPLIFY = FALSE)
}

rm(z)

######################################################################################

## call binding events

bs_dir <- file.path(base_dir,"binding")
check_create(bs_dir)

bs_dirs <- lapply(what,function(x)file.path(bs_dir,x))
lapply(bs_dirs,check_create)

bs_dirs <- file.path(bs_dirs,paste0("FDR",fdr * 100))
lapply(bs_dirs,check_create)

bs_dirs <- file.path(bs_dirs,paste0("G_",maxComp))
lapply(bs_dirs,check_create)


call_sites <- function(peak_file,read_file,out_file,fl,g,pet,mc)
{
  dp <- dpeakRead(peakfile = peak_file,readfile = read_file, fileFormat = "bam",
    PET = pet, fragLen = fl,parallel = TRUE, nCore = mc)
  dp <- dpeakFit(dp,nCore = mc)
  export(dp,type = "bed",filename = out_file)
  message("*********************")
}

dpeak_sites_wrap <- function(peak_dir,read_dir,out_dir,dt,what,fl,g,mc)
{
  if(what == "pet"){
    pet = TRUE
  }else{
    pet = FALSE
  }

  peaks <- list.files(peak_dir)
  reads <- sample_read_files(read_dir,dt,pet)

  if(is.unsorted(peaks))peaks <- sort(peaks)
  if(is.unsorted(reads))reads <- sort(reads)

  outs <- sapply(strsplit(reads,".",fixed = TRUE),function(x)x[1])
  outs <- file.path(out_dir,paste0(outs,".bed"))
                        
  peaks <- file.path(peak_dir,peaks)
  reads <- file.path(read_dir,reads)
  
  sites <- mapply(call_sites,
    peaks,reads,outs,MoreArgs = list(fl,g,pet,mc),SIMPLIFY = FALSE)

  u <- gc()
}

if(flag_binding){
  z <- mapply(dpeak_sites_wrap,peak_dirs,chip_dirs,bs_dirs,char,what,
    MoreArgs = list(fl,maxComp,mc),SIMPLIFY = FALSE)
}

rm(z)

######################################################################################
