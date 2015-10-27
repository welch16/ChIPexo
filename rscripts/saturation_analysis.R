
rm(list = ls())

library(mosaics)
library(data.table)
library(GenomicAlignments)
library(dpeak)
## library(devtools)
## load_all("~/Desktop/Docs/Code/dpeak")


######################################################################################

## Condition table

source("R/base_edsn.R")

exo_char <- edsn_tab("exo")
pet_char <- edsn_tab("pet")
set_char <- edsn_tab("set")

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif20min"
bs <- 150
fl <- 150
fdr <- .001
thresh <- 10
mc <- 16
maxComp <- 5
seedset <- c( 23456, 34567, 45678, 56789,45987 )

### we have 5 different seeds
k <- 2
flag_sample <- FALSE
flag_bins <- FALSE
flag_peaks <- FALSE
flag_binding <- TRUE

## notes on seeds:
##
## 23456 works on peaks
##
## 12345 dont 

exo <- exo_char[ip == tf & condition == rif ]
pet <- pet_char[ip == tf & condition == rif ]
set <- set_char[ip == tf & condition == rif ]

n_val <- seq(5e4,2e5,by = 1e4)

base_dir <- "/p/keles/ChIPexo/volume6/saturation"
folder <- paste(tf,rif,sep = "_")
what <- c("exo","pet","set")

base_dir <- file.path(base_dir,folder)
if(!dir.exists(base_dir))dir.create(base_dir)

base_dir <- file.path(base_dir,paste0("seed",k))
if(!dir.exists(base_dir))dir.create(base_dir)

##################################h####################################################

## Sample

dir1 <- "/p/keles/ChIPexo/volume3/LandickData"
exo_reads_dir <- file.path(dir1,"ChIPexo")
pet_reads_dir <- file.path(dir1,"ChIPseq_PET")
set_reads_dir <- file.path(dir1,"ChIPseq_SET")

get_reads_files <- function( dir,dt,pet)
{
  edsn <- dt[,(edsn)]
  files <- list.files(dir)
  files <- files[sapply(edsn,function(x)grep(x,files))]
  if(pet){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

exo_reads_files <- get_reads_files(exo_reads_dir,exo,FALSE)
pet_reads_files <- get_reads_files(pet_reads_dir,pet,TRUE)
set_reads_files <- get_reads_files(set_reads_dir,set,FALSE)

filter_factory <- function(want){
  list(KeepQname = function(x) x$qname %in% want)
}

sample_dir <- file.path(base_dir,"sample")
if(!dir.exists(sample_dir))dir.create(sample_dir)

sample_reads <- function(reads_file,in_dir,out_dir,n_val,pet)
{
  ff <- file.path(in_dir,reads_file)
  edsn <- strsplit(reads_file,"_")[[1]][1]

  ff_out <- file.path(out_dir,
    sapply(n_val,function(x)paste0(edsn,"_sample",as.character(x/1000),"K.bam")))
  
  if(pet){
    param <- ScanBamParam(what = c("qname"))
    pairs <- readGAlignmentPairs(ff , param = param)
    qn <- mcols(left(pairs))$qname  ## since left's and right's qnames are identicall
    want_list <- lapply(n_val,function(x)sample(qn,x))
  }else{
    param <- ScanBamParam(what = "qname")
    reads <- readGAlignments(ff , param = param)
    qn <- mcols(reads)$qname
    want_list <- lapply(n_val,function(x)sample(qn,x))
  }
    
  out <- mapply(function(want,dest,pp){
    filter <- FilterRules(filter_factory(want))
    ffo <- filterBam(ff,dest,filter = filter, param = pp)    
  },want_list,ff_out,MoreArgs = list(param),SIMPLIFY = FALSE)

  return(out)
}

sample_dirs <- sapply(what,function(x)file.path(sample_dir,x))
lapply(sample_dirs,function(x)if(!dir.exists(x))dir.create(x))

if(flag_sample){
  set.seed(seedset[k])
  lapply(exo_reads_files,sample_reads,exo_reads_dir,sample_dirs[[1]],n_val,FALSE)

  set.seed(seedset[k])
  lapply(pet_reads_files,sample_reads,pet_reads_dir,sample_dirs[[2]],n_val,TRUE)

  set.seed(seedset[k])
  lapply(set_reads_files,sample_reads,set_reads_dir,sample_dirs[[3]],n_val,FALSE)

}

######################################################################################

### build_bins

bin_dir <- file.path(base_dir,"bins")
if(!dir.exists(bin_dir))dir.create(bin_dir)

bin_dirs <- lapply(what,function(x)file.path(bin_dir,x))
lapply(bin_dirs,function(x)if(!dir.exists(x))dir.create(x))

sample_read_files <- function(in_dir,dt)
{
  files <- list.files(in_dir)
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

create_bins <- function(in_dir,out_dir,dt,pet,bs,fl)
{
  files <- sample_read_files(in_dir,dt)
  lapply(file.path(in_dir,files),constructBins,
         fileFormat = "bam",
         outfileLoc = out_dir,
         PET = pet,
         fragLen = fl,
         binSize = bs)                  
}

if(flag_bins){
  mapply(create_bins,sample_dirs,bin_dirs,list(exo,pet,set),list(FALSE,TRUE,FALSE),
    MoreArgs = list(bs,fl),SIMPLIFY = FALSE)
}

######################################################################################

## call peaks

peak_dir <- file.path(base_dir,"peaks")
if(!dir.exists(peak_dir))dir.create(peak_dir)

peak_dirs <- lapply(what,function(x)file.path(peak_dir,x))
lapply(peak_dirs,function(x)if(!dir.exists(x))dir.create(x))

peak_dirs <- file.path(peak_dirs,paste0("FDR",fdr*100))
lapply(peak_dirs,function(x)if(!dir.exists(x))dir.create(x))

extra <- list()
extra[["exo"]] <- file.path("/p/keles/genome_data/EColi_U00096.2",
  c("mappability","GC","N"),"bin",
  c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
    "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))
extra[["pet"]] <- "/p/keles/ChIPexo/volume6/inputs/pet/edsn1369_042814_qc.filter.bam_bin150.txt"
extra[["set"]] <- "/p/keles/ChIPexo/volume6/inputs/set/edsn1369_042814_qc.filter.bam_fragL150_bin150.txt"

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
if(!dir.exists(bs_dir))dir.create(bs_dir)

bs_dirs <- lapply(what,function(x)file.path(bs_dir,x))
lapply(bs_dirs,function(x)if(!dir.exists(x))dir.create(x))

bs_dirs <- file.path(bs_dirs,paste0("G_",maxComp))
lapply(bs_dirs,function(x)if(!dir.exists(x))dir.create(x))


call_sites <- function(peak_file,read_file,out_file,fl,g,pet,mc)
{
  dp <- dpeakRead(peakfile = peak_file,readfile = read_file, fileFormat = "bam",
    PET = pet, fragLen = fl,parallel = TRUE, nCore = mc)
  dp <- dpeakFit(dp,nCore = mc,maxComp = g)
  export(dp,type = "bed",filename = out_file)
}

dpeak_sites_wrap <- function(peak_dir,read_dir,out_dir,what,fl,g,mc)
{
  if(what == "pet"){
    pet = TRUE
  }else{
    pet = FALSE
  }

  peaks <- list.files(peak_dir)
  reads <- list.files(read_dir)
  reads <- reads[grep("bai",reads,invert =TRUE)]

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

##dpeak_sites_wrap(peak_dirs[2],sample_dirs[2],bs_dirs[2],"pet",fl,maxComp,mc)

if(flag_binding){
  z <- mapply(dpeak_sites_wrap,peak_dirs,sample_dirs,bs_dirs,c("exo","pet","set"),
    MoreArgs = list(fl,maxComp,mc),SIMPLIFY = FALSE)
}

rm(z)

######################################################################################
