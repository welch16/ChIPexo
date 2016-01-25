#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    call_peaks.R - Using mosaics call peaks for all the bin files in the indir.

  Arguments:


   -- indir

      The directory where the bam files with the sampled reads are located.

   -- outdir

      The location where the bins are going to be saved.

   -- fdr

      Value use to filter the peaks, considered as a percent (5% as 5)
 
   -- what

      Character value of the type of additional data used for the model. Can be: exo, pet or set.
      
   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 4)

library(GenomicAlignments)
library(parallel)
library(mosaics)
library(data.table)

indir <- args[1]
outdir <- args[2]
fdr <- as.numeric(args[3]) / 100
what <- args[4]

indir <- file.path(indir,"bins")
nCore <- detectCores()
files <- list.files(indir)

if(what == "exo"){
  extra <- file.path("/p/keles/genome_data/EColi_U00096.2",
    c("mappability","GC","N"),
    "bin",
    c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))
}else if(what == "pet"){
  extra <- "/p/keles/ChIPexo/volume6/resolution/inputs/ChIPseq_PET/edsn1369_Input.sort.bam_bin150.txt"
}else{
  extra <- "/p/keles/ChIPexo/volume6/resolution/inputs/ChIPseq_SET/edsn1369_Input.sort.bam_fragL150_bin150.txt"
}

call_peaks <- function(binfile,extra,what,fdr)
{
  if(what == "exo"){
    type <- c("chip","M","GC","N")
  }else{
    type <- c("chip","input")
  }

  bins <- readBins(type = type, fileName = c(binfile,extra))
  fit <- try(mosaicsFit(bins),silent = TRUE)
  if(class(fit) == "try-error"){
    peaks <- data.table("nopeaks")
  }else{
    mod <- ifelse(fit@bic1S < fit@bic2S ,"1S","2S")
    mp <- mosaicsPeak(fit,signalModel = mod,thres = 100)
    peaks <- data.table(print(mp))
  }
  return(peaks)
}

peaks <- mclapply(file.path(indir,files),
  call_peaks,extra,what,fdr,mc.cores = nCore)

outfiles <- strsplit(files,".",fixed = TRUE)
outfiles <- sapply(outfiles,function(x)x[1])
outfiles <- paste0(outfiles,"_peaks.txt")

a <- mapply(write.table,
  peaks,file = file.path(outdir,outfiles),
  MoreArgs = list(row.names = FALSE,col.names = FALSE,quote = FALSE))            




