#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    call_binding_sites.R - Deconvolve the peak signal into binding events.

  Arguments:


   -- readsfile

      The bam file where the reads are stored.

   -- peaksfile

      A fixed peak file that is gonna be used to find binding events in its peaks.

   -- outdir

      Directory where the output is gonna be saved.

   -- isPET

      Logical variable indicating if the reads are Paired Ended or not.

   -- maxG

      Integer value indicating what is the maximum number of binding events per region.
      
   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}


stopifnot(length(args) == 5)


library(GenomicAlignments)
library(parallel)
library(dpeak)
library(data.table)

readsfile <- args[1]
peaksfile <- args[2]
outdir <- args[3]
isPET <- as.logical(args[4])
maxG <- as.numeric(args[5])

idx <- strsplit(basename(readsfile),"_",fixed = TRUE)[[1]]

fline <- scan(peaksfile,what = character(),nlines = 1)

outfile <- basename(gsub(".bam","_sites.txt",readsfile))
outfile <- file.path(outdir,outfile)
                    

if(length(fline) == 1 & fline == "nopeaks"){
  write.table("nobinding",file = outfile,col.names = FALSE,row.names = FALSE,quote = FALSE)
}else{
  nCore <- detectCores()
  dpeak <- dpeakRead(peakfile = peaksfile,readfile = readsfile,fileFormat ="bam",PET = isPET,nCore = nCore)
  fit <- dpeakFit(dpeak, maxComp = maxG,nCore = nCore)
   tt <- tempfile(fileext = "bed")
  export(fit,type = ".bed",filename = tt)
  dt <- data.table(read.table(tt,skip = 1))
  setnames(dt,names(dt),c("chrID","start","end","peakID","strength"))
  write.table(dt,file = outfile,quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
}
