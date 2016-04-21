#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    SCC_curve.R - Calculates the SCC curve for a given ChIP (seq,exo
      or nexus) sample

  Arguments:

   -- bamfile

      File in bam format with the aligned reads used to calculate the curve

   -- outfile

      Name of the file, where the output is gonna be saved

   -- sizefile

      File without header and with two column indicating the
      respective chromosome and it's length. If the data was aligned to
      any of the dm3, hg19, mm9 or mm10 genomes, then it loads the file
      automatically by using dm3, hg19, etc.

   -- isPET

      Boolean variable indicating if the reads in the bamfile are paired

   -- maxShift

      Max possible shift for the curve, i.e. the output is gonna be a
      data.table with format [shift = 1:maxShift , cross.corr]

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 5)

bamfile <- args[1]
outfile <- args[2]
sizefile <- args[3]
isPET <- as.logical(args[4])
maxShift <- as.numeric(args[5])

stopifnot(file.exists(bamfile))
stopifnot(maxShift > 0)

library(parallel)
library(data.table)
library(GenomicAlignments)
library(devtools)

load_all("~/Desktop/Docs/Code/ChIPUtils")

mc <- detectCores()

reads <- create_reads(bamfile,isPET)

if(tolower(sizefile) %in% c("hg19","mm9","mm10","dm3")){
  sizedir <- system.file("extdata","chrom.sizes", package = "ChIPUtils")
  sizefiles <- list.files(sizedir)
  sizefile <- sizefiles[grep(sizefile,sizefiles)]
  sizefile <- file.path(sizedir,sizefile)
  rm(sizedir,sizefiles)
}

sizes <- data.table(read.table(sizefile,header = FALSE))

scc <- strand_cross_corr(reads,shift = 1:maxShift,
   chrom.sizes = sizes,parallel = TRUE)

write.table(format(scc,digits = 6),file = outfile,quote = FALSE,
   sep = "\t",row.names = FALSE,col.names = TRUE)            


