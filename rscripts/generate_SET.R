#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    generate_SET.R - Samples a SET file out of a PET file a write it in the same location as file

  Arguments:

    -- file

      A bam file, with the reads of PET ChIP-Seq experiment

   -- help

      Show the help file

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

library(GenomicAlignments)

file <- args[1]
dir <- dirname(file)
out_file <- gsub(".sort.bam",".SET.bam",file)

filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}

param1 <- ScanBamParam(scanBamFlag(isFirstMateRead = TRUE, isSecondMateRead = FALSE))
param2 <- ScanBamParam(scanBamFlag(isFirstMateRead = FALSE, isSecondMateRead = TRUE))

fmate <- readGAlignments(file, param = param1,use.names = TRUE)
nms <- names(fmate)

set.seed(123)
nn <- rbinom(1,prob = .5,size = length(nms))

idx <- sample(1:length(nms),nn)

filter1 <- FilterRules(filter_factory(nms[idx]))
filter2 <- FilterRules(filter_factory(nms[-idx]))

f1 <- tempfile("reads",fileext=".bam")
f2 <- tempfile("reads",fileext=".bam")

param1 <- ScanBamParam(what = "qname",scanBamFlag(isFirstMateRead = TRUE, isSecondMateRead = FALSE))
param2 <- ScanBamParam(what = "qname",scanBamFlag(isFirstMateRead = FALSE, isSecondMateRead = TRUE))

ff1 <- filterBam(file,f1,filter = filter1,param = param1)
ff2 <- filterBam(file,f2,filter = filter2,param = param2)

ff <- mergeBam(c(ff1,ff2),out_file)


