
rm(list = ls())

library(mosaics)
library(data.table)
library(GenomicAlignments)
library(dpeak)

######################################################################################

## Condition table

edsn_tab <- function(what){
  stopifnot(what %in% c("exo","pet","set"))
  if(what == "exo"){
    edsn <- as.character(931:938)
    ip <- rep(c("Sig70","SigmaS"),4)
    condition <- rep(c("exp","stat"),each = 4)
    repl <- rep( rep(1:2,each =2),2)
    dt1 <- data.table(edsn,ip,condition,repl)
    edsn <- as.character(1310:1321)
    ip <- rep(c("Beta","Sig70","BetaPF"),4)
    condition <- rep( rep(c("rif0min","rif20min"),each = 3),2)
    repl <- rep(rep(1:2,each = 6))
    dt2 <- data.table(edsn,ip,condition,repl)
    dt <- rbind(dt1,dt2)    
  }else{
    edsn <- as.character(1396:1403)
    ip <- rep(c("Sig70","BetaPF"),4)
    condition <- rep(c("rif0min","rif20min"),each = 2, 2)
    repl <- rep(1:2,each = 4)
    dt <- data.table(edsn,ip,condition,repl)      
  }
  return(dt)
}
  
exo_char <- edsn_tab("exo")
pet_char <- edsn_tab("pet")
set_char <- edsn_tab("set")

#############1#########################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif0min"
bs <- 150
fl <- 150
fdr <- .01
thresh <- 10
mc <- 16
maxComp <- 1

flag_bins <- FALSE
flag_peaks <- TRUE
flag_binding <- TRUE
flag_input <- FALSE

exo <- exo_char[ip == tf & condition == rif ]
pet <- pet_char[ip == tf & condition == rif ]
set <- set_char[ip == tf & condition == rif ]

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

base_dir <- "/p/keles/ChIPexo/volume6/condition"
folder <- paste(tf,rif,sep = "_")
what <- c("exo","pet","set")

base_dir <- file.path(base_dir,folder)
check_create(base_dir)

##################################h####################################################

## Input data bins for ChIPseq (PET and SET)

inputs_dir <- "/p/keles/ChIPexo/volume6/inputs"
inputs_dirs <- lapply(what[-1],function(x)file.path(inputs_dir,x))
  
if(flag_input){
  lapply(inputs_dirs,check_create)
  constructBins(file.path(pet_dir,"edsn1369_042814_qc.filter.bam"),
    fileFormat = "bam",outfileLoc = inputs_dirs[[1]],
    PET = TRUE,fragLen = fl,binSize = bs)
  constructBins(file.path(set_dir,"edsn1369_042814_qc.filter.bam"),
    fileFormat = "bam",outfileLoc = inputs_dirs[[2]],
    PET = FALSE,fragLen = fl,binSize = bs)                              
}

######################################################################################

### build_bins

chip_dirs <- list()
chip_dirs[["exo"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
chip_dirs[["pet"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
chip_dirs[["set"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

bin_dir <- file.path(base_dir,"bins")
check_create(bin_dir)

bin_dirs <- lapply(what,function(x)file.path(bin_dir,x))
lapply(bin_dirs,check_create)

sample_read_files <- function(in_dir,dt,pet)
{
  files <- list.files(in_dir)
  edsn <- dt[,(edsn)]
  files <- files[sapply(edsn,function(x)grep(x,files))]
  if(pet){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

create_bins <- function(in_dir,out_dir,dt,pet,bs,fl)
{
  files <- sample_read_files(in_dir,dt,pet)
  lapply(file.path(in_dir,files),constructBins,
         fileFormat = "bam",
         outfileLoc = out_dir,
         PET = pet,
         fragLen = fl,
         binSize = bs)                  
}

if(flag_bins){
  mapply(create_bins,chip_dirs,bin_dirs,list(exo,pet,set),list(FALSE,TRUE,FALSE),
    MoreArgs = list(bs,fl),SIMPLIFY = FALSE)
}

######################################################################################

## call peaks

peak_dir <- file.path(base_dir,"peaks")
check_create(peak_dir)

peak_dirs <- lapply(what,function(x)file.path(peak_dir,x))
lapply(peak_dirs,check_create)

peak_dirs <- file.path(peak_dirs,paste0("FDR",fdr*100))
lapply(peak_dirs,check_create)

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
check_create(bs_dir)

bs_dirs <- lapply(what,function(x)file.path(bs_dir,x))
lapply(bs_dirs,check_create)

bs_dirs <- file.path(bs_dirs,paste0("FDR",fdr * 100))
lapply(bs_dirs,check_create)

bs_dirs <- file.path(bs_dirs,paste0("G_",maxComp))
lapply(bs_dirs,check_create)


call_sites <- function(peak_file,read_file,out_file,fl,g,pet,mc)
{
  dp <- dpeakRead(peakfile = peak_file,readfile = read_file, fileFormat = "bam",
    PET = pet, fragLen = fl,parallel = TRUE, nCore = mc)
  dp <- dpeakFit(dp,nCore = mc)
  export(dp,type = "bed",filename = out_file)
}

dpeak_sites_wrap <- function(peak_dir,read_dir,out_dir,dt,what,fl,g,mc)
{
  if(what == "pet"){
    pet = TRUE
  }else{
    pet = FALSE
  }

  peaks <- list.files(peak_dir)
  reads <- sample_read_files(read_dir,dt,pet)

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

if(flag_binding){
  z <- mapply(dpeak_sites_wrap,peak_dirs,chip_dirs,bs_dirs,list(exo,pet,set),c("exo","pet","set"),
    MoreArgs = list(fl,maxComp,mc),SIMPLIFY = FALSE)
}

rm(z)

######################################################################################
