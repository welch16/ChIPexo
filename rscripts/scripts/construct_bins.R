#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    construct_bins.R - Samples a SET file out of a PET file a write it in the same location as file.

  Arguments:

   -- indir

      The directory where the bam files with the sampled reads are located.

   -- outdir

      The location where the bins are going to be saved.

   -- fragLen

      The fragment length used to extend the fragments in case of SET experiment. This parameter
      is ignored in the PET case.

   -- binSize

      The bin size used to construct the bins.

   -- isPET

      Logical value indicating if the reads are Paired End or Single End.
      
   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 5)

library(GenomicAlignments)
library(parallel)
library(mosaics)

indir <- args[1]
outdir <- args[2]
fragLen <- as.numeric(args[3])
binSize <- as.numeric(args[4])
isPET <- as.logical(args[5])

nCore <- detectCores()
files <- list.files(indir)

files <- files[grep("bai",files,invert = TRUE)]

A <- mclapply(file.path(indir,files),
  constructBins,
  fileFormat = "bam",
  outfileLoc = outdir,
  PET = isPET,fragLen = fragLen,binSize = binSize,mc.cores = nCore)           
