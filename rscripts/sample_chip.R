
rm(list = ls())

library(GenomicAlignments)
library(Rsamtools)


base_dir <- "/p/keles/ChIPexo/volume3"
input_dirs <- list()
input_dirs[["landick_exo"]] <- file.path(base_dir,"LandickData/ChIPexo/")
input_dirs[["landick_pet"]] <- file.path(base_dir,"LandickData/ChIPseq_PET")

get_files <- function(dr)
{
  files <- list.files(dr)
  files <- files[grep("sort",files)]
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

files <- lapply(input_dirs,get_files)
nValues <- floor( seq(1e4,1e6,length.out = 25))

filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}

sample_exo <- function(file,in_dir,n_val,out_dir)
{

  ff <- file.path(in_dir,file)
  message(ff)

  ff_out <- file.path(out_dir,
       sapply(1:length(n_val),function(x) gsub(".sorted.bam",paste0(".sample",x ,".bam"),file)))

  param <- ScanBamParam(what = "qname")
  reads <- readGAlignments(ff , param = param)

  qn <- mcols(reads)$qname
  want_list <- lapply(n_val,function(x)sample(qn,x))

  out <- mapply(function(want,dest,pp){
    filter <- FilterRules(filter_factory(want))
    ffo <- filterBam(ff,dest,filter = filter, param = pp)
    
  },want_list,ff_out,MoreArgs = list(param),SIMPLIFY = FALSE)

  
}


sample_pet <- function(file,in_dir,n_val,out_dir)
{

  ff <- file.path(in_dir,file)
  message(ff)


  ff_out <- file.path(out_dir,
       sapply(1:length(n_val),
         function(x) gsub(".sorted.bam",paste0(".sample",x ,".bam"),file)))

  param <- ScanBamParam(what = c("qname"))
  pairs <- readGAlignmentPairs(ff , param = param)

  qn <- mcols(left(pairs))$qname  ## since left's and right's qnames are identicall
  want_list <- lapply(2*n_val,function(x)sample(qn,x))

  out <- mapply(function(want,dest,pp){
    filter <- FilterRules(filter_factory(want))
    ffo <- filterBam(ff,dest,filter = filter, param = pp)    
  },want_list,ff_out[[1]],MoreArgs = list(param),SIMPLIFY = FALSE)
  
}


sample_set <- function(file,in_dir,n_val,out_dir)
{
  ff <- file.path(in_dir,file)
  message(ff)

  ff_out <- file.path(out_dir,
     sapply(1:length(n_val),
         function(x) gsub(".sorted.bam",paste0(".sample",x ,".bam"),file)))

  mate1 <- sapply(1:length(n_val),function(x) tempfile("reads",fileext= ".bam"))  
  mate2 <- sapply(1:length(n_val),function(x) tempfile("reads",fileext= ".bam"))

  param1 <- ScanBamParam(what = c("qname"),scanBamFlag(isFirstMateRead = TRUE, isSecondMateRead = FALSE))
  fmate <- readGAlignments(ff , param = param1)
  qn1 <- mcols(fmate)[["qname"]]

  param2 <- ScanBamParam(what = c("qname"),scanBamFlag(isFirstMateRead = FALSE, isSecondMateRead = TRUE))
  smate <- readGAlignments(ff , param = param2)
  qn2 <- mcols(smate)[["qname"]]

  want_list <- lapply(n_val,function(x){
    nn <- rbinom( 1 , prob = .5 , size = x)
    out <- list()
    out[["mate1"]] <- sample(qn1,nn)
    out[["mate2"]] <- sample(qn2,x - nn)
    return(out)
  })

  out1 <- mapply(function(want,dest,pp){
    filter <- FilterRules(filter_factory(want[["mate1"]]))
    ffo <- filterBam(ff,dest,filter = filter, param = pp)    
  },want_list,mate1,MoreArgs = list(param1),SIMPLIFY = FALSE)

  out2 <- mapply(function(want,dest,pp){
    filter <- FilterRules(filter_factory(want[["mate2"]]))
    ffo <- filterBam(ff,dest,filter = filter, param = pp)    
  },want_list,mate2,MoreArgs = list(param2),SIMPLIFY = FALSE)

 out <- mapply(function(m1,m2,out){
   files <- c(m1,m2)
   mergeBam(files, out,overwrite = TRUE)
 },mate1,mate2,ff_out,SIMPLIFY = FALSE)

  lapply(out,indexBam)
  
  
}

out_dir <- "/p/keles/ChIPexo/volume6/saturation/Landick/ChIPexo"
exo <- lapply(files[[1]],sample_exo,input_dirs[[1]],nValues,out_dir)

seq_dirs <- file.path("/p/keles/ChIPexo/volume6/saturation/Landick",c("ChIPseq_PET","ChIPseq_SET"))
pet <- lapply(files[[2]],sample_pet,input_dirs[[2]],nValues,seq_dirs[[1]])

set <- lapply(files[[2]],sample_set,input_dirs[[2]],nValues,seq_dirs[[2]])



 
