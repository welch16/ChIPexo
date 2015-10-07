
rm(list = ls())

library(GenomicAlignments)
library(Rsamtools)


base_dir <- "/p/keles/ChIPexo/volume3"
in_dir <- file.path(base_dir,"LandickData/ChIPseq_PET")
out_dir <- file.path(base_dir,"LandickData/ChIPseq_SET")

get_files <- function(dr,PET = FALSE)
{
  files <- list.files(dr)
  files <- files[grep(ifelse(PET,"filter","sort"),files)]
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

files <- get_files(in_dir,TRUE)

filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}

sample_set <- function(file,in_dir,out_dir)
{
  ff <- file.path(in_dir,file)
  message(ff)

  ff_out <- file.path(out_dir,file)

  mate1 <- tempfile("reads",fileext= ".bam")
  mate2 <- tempfile("reads",fileext= ".bam")

  param <- ScanBamParam(what = c("qname"))
  pairs <- readGAlignmentPairs(ff, param = param)

  N <- length(left(pairs))
  
  M1 <- rbinom(1,prob = .5,size = N)
  q1 <- mcols(left(pairs))[["qname"]]
  idx <- sample(N, M1)

  want1 <- q1[idx]
  want2 <- q1[-idx]

  filter1 <- FilterRules(filter_factory(want1))
  filter2 <- FilterRules(filter_factory(want2))
  
  param1 <- ScanBamParam(what = c("qname"),scanBamFlag(isFirstMateRead = TRUE, isSecondMateRead = FALSE))
  param2 <- ScanBamParam(what = c("qname"),scanBamFlag(isFirstMateRead = FALSE, isSecondMateRead = TRUE))

  out1 <- filterBam(ff,mate1,filter = filter1,param = param1)
  out2 <- filterBam(ff,mate2,filter = filter2,param = param2)

  out <- mergeBam(c(mate1,mate2),ff_out,overwrite = TRUE)

  indexBam(out)

  
}

set.seed(1234321)

set <- lapply(files,sample_set,in_dir,out_dir)



 
