#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    generate_SET.R - Samples a SET file out of a PET file a write it in the same location as file

  Arguments:

    -- infile

      The bam file where the reads are stored. It is assumed that the file is named
      with the suffix .sort.bam

   -- outdir

      The directory where the sampled bam files are gonna be saved.

   -- minN, maxN and inc

      The minimun, maximum and the increment used to sample the reads. If the reads are paired ended,
      we are going to sample  2 * N reads instead of N to complete the pairs.

   -- isPET

      Logical value indicating if the reads are Paired End or Single End.

   -- seed

      Seed used to sample the reads.

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 7)

infile <- args[1]
outdir <- args[2]
minN <- as.numeric(args[3])
maxN <-as.numeric(args[4])
inc <- as.numeric(args[5])
isPET <- as.logical(args[6])
seed <- as.numeric(args[7])

library(GenomicAlignments)
library(parallel)

nCore <- detectCores()

set.seed(seed)

samps <- seq(minN,maxN,by = inc)

filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}


sample_bam <- function(outfile,N,file,isPET){

  param <- NULL
  if(isPET){
    pairs <- readGAlignmentPairs(file,param = NULL,use.names = TRUE)
    want <- sample(names(pairs),N)
  }else{
    reads <- readGAlignments(file,param = NULL,use.names = TRUE)
    want <- sample(names(reads),N)
  }

  filter <- FilterRules(filter_factory(want))
  out <- filterBam(file,outfile,filter = filter,param = ScanBamParam(what = "qname"))
}

outfiles <- gsub(".sort.bam","",basename(infile))
outfiles <- paste0(outfiles,"_samp",samps,"_",seed,".bam")
outfiles <- file.path(outdir,outfiles)

a = mcmapply(sample_bam,outfiles,samps,MoreArgs = list(infile,isPET),
  SIMPLIFY = FALSE,mc.set.seed = FALSE,mc.cores = nCore)





