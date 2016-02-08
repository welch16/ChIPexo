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

   -- outfile

      The directory where the sampled bam files are gonna be saved.

   -- N

      Number of sampled reads. (For PET files is gonna sample N / 2)


   -- isPET

      Logical value indicating if the reads are Paired End or Single End.

   -- seed

      Seed used to sample the reads.

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 5)

infile <- args[1]
outfile <- args[2]
N <- as.numeric(args[3])
isPET <- as.logical(args[4])
seed <- as.numeric(args[5])

library(GenomicAlignments)

set.seed(seed)


filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}


sample_bam <- function(outfile,N,file,isPET){

  param <- NULL
  if(isPET){
    pairs <- readGAlignmentPairs(infile,param = NULL,use.names = TRUE)
    want <- sample(names(pairs),N)
  }else{
    reads <- readGAlignments(file,param = NULL,use.names = TRUE)
    want <- sample(names(reads),N)
  }

  filter <- FilterRules(filter_factory(want))
  out <- filterBam(file,outfile,filter = filter,param = ScanBamParam(what = "qname"))
}

if(isPET){
  N <- floor(N/2)
}

out <- sample_bam(outfile,N,infile,isPET)
