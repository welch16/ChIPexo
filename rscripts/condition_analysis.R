
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

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif20min"
bs <- 150
fl <- 150
fdr <- .25
thresh <- 10
mc <- detectCores()
g <- 5

exo <- exo_char[ip == tf & condition == rif ]
pet <- pet_char[ip == tf & condition == rif ]
set <- set_char[ip == tf & condition == rif ]

base_dir <- "/p/keles/ChIPexo/volume6/condition"
folder <- paste(tf,rif,sep = "_")

base_dir <- file.path(base_dir,folder)

if(!dir.exists(base_dir)){
  dir.create(base_dir)
}

what <- c("exo","pet","set")
bases <- lapply(what,function(x)file.path(base_dir,x))
lapply(bases,function(x)if(!dir.exists(x))dir.create(x))

######################################################################################

### build_bins
exo_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
pet_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
set_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

read_files <- function(in_dir,dt,pet)
{
  files <- list.files(in_dir)
  edsn <- dt[,(edsn)]    
  files <- sapply(edsn,function(x,files)files[grep(x,files)],files)
  if(pet){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}


create_bins <- function(in_dir,out_dir,dt,bs,fl,pet)
{
  files <- read_files(in_dir,dt,pet)
  lapply(file.path(in_dir,files),constructBins,
         fileFormat = "bam",
         outfileLoc = out_dir,
         PET = pet,
         fragLen = fl,
         binSize = bs)                  
}

create_bins(exo_dir,bases[[1]],exo,bs,fl,FALSE)
create_bins(pet_dir,bases[[2]],pet,bs,fl,TRUE)
create_bins(set_dir,bases[[3]],set,bs,fl,FALSE)

######################################################################################

## call peaks

### chip exo

exo_bins <- list.files(bases[[1]])
input_dir <- "/p/keles/genome_data/EColi_U00096.2/"

read_bins_wrap <- function(file,map_file,gc_file,n_file,opt)
{
  readBins(opt,c(file, map_file,gc_file,n_file))
}

mosaics_peak_wrap <- function(fit,FDR,binsize,thres)
{
  if ( fit@bic2S <= fit@bic1S ) {
    peak <- mosaicsPeak( fit, signalModel="2S", FDR=FDR, maxgap=binsize, thres=thres )
  } else {
    peak <- mosaicsPeak( fit, signalModel="1S", FDR=FDR, maxgap=binsize, thres=thres )
  }
  return(peak)
}

exo_bins <- lapply(file.path(bases[[1]],exo_bins),
  read_bins_wrap,
  map_file = file.path(input_dir,
    "mappability/bin/fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt"),
  gc_file = file.path(input_dir,
    "GC/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt"),
  n_file = file.path(input_dir,
    "N/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"),
  opt = c("chip","M","GC","N"))
                             
exo_fits <- mclapply(exo_bins,mosaicsFit,mc.cores = 2)
exo_peaks <- mclapply(exo_fits,mosaics_peak_wrap,
                        fdr,bs,thresh,mc.cores = 2)
exo_peaks <- lapply(exo_peaks,function(x)data.table(x@peakList))

exo_peaks_loc <- file.path(bases[[1]],paste0("exo_peaks_",1:2,".txt"))
mapply(write.table,
       exo_peaks,exo_peaks_loc,
       MoreArgs = list(col.names =FALSE,row.names = FALSE,quote = FALSE))

### chip seq pet

pet_bins <- list.files(bases[[2]])

read_bins_wrap <- function(file,input_file,opt)
{
  readBins(opt,c(file, input_file))
}


pet_bins <- lapply(file.path(bases[[2]],pet_bins),
  read_bins_wrap,
  input_file = "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_PET/bins/edsn1369_042814_qc.filter.bam_bin150.txt",
  opt = c("chip","input"))
                             
pet_fits <- mclapply(pet_bins,mosaicsFit,mc.cores = 2)
pet_peaks <- mclapply(pet_fits,mosaics_peak_wrap,
                        fdr,bs,thresh,mc.cores = 2)
pet_peaks <- lapply(pet_peaks,function(x)data.table(x@peakList))

pet_peaks_loc <- file.path(bases[[2]],paste0("pet_peaks_",1:2,".txt"))
mapply(write.table,
       pet_peaks,pet_peaks_loc,
       MoreArgs = list(col.names =FALSE,row.names = FALSE,quote = FALSE))

### chip seq set

set_bins <- list.files(bases[[3]])

set_bins <- lapply(file.path(bases[[3]],set_bins),
  read_bins_wrap,
  input_file = "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_SET/bins/edsn1369_042814_qc.filter.bam_fragL150_bin150.txt",
  opt = c("chip","input"))
                             
set_fits <- mclapply(set_bins,mosaicsFit,mc.cores = 2)
set_peaks <- mclapply(set_fits,mosaics_peak_wrap,
                        fdr,bs,thresh,mc.cores = 2)
set_peaks <- lapply(set_peaks,function(x)data.table(x@peakList))

set_peaks_loc <- file.path(bases[[3]],paste0("set_peaks_",1:2,".txt"))
mapply(write.table,
       set_peaks,set_peaks_loc,
       MoreArgs = list(col.names =FALSE,row.names = FALSE,quote = FALSE))

rm(exo_bins,pet_bins,set_bins,exo_fits,pet_fits,set_fits,exo_peaks,set_peaks,pet_peaks,
   opt,input_dir,what,folder,base_dir)

######################################################################################

## call binding events

### chip exo binding sites

exo_read_loc <- file.path(exo_dir,read_files(exo_dir,exo,FALSE))

exo_dpeak <- mapply(dpeakRead,exo_peaks_loc,exo_read_loc,
  MoreArgs = list(fileFormat = "bam",PET = FALSE,fragLen = fl,parallel = TRUE,nCore = mc),SIMPLIFY =FALSE)

exo_dpeak <- lapply(exo_dpeak,dpeakFit,maxComp = g,nCore = mc)

exo_bs_loc <- file.path(bases[[1]],paste0("exo_bs_g",g,"_",1:2,".bed"))

for(i in 1:2){
  export(exo_dpeak[[i]],type = "bed",filename = exo_bs_loc[i])
}

## chip set pet binding sites

pet_read_loc <- file.path(pet_dir,read_files(pet_dir,pet,TRUE))

pet_dpeak <- mapply(dpeakRead,pet_peaks_loc,pet_read_loc,
  MoreArgs = list(fileFormat = "bam",PET = TRUE,fragLen = fl,parallel = TRUE,nCore = mc),SIMPLIFY =FALSE)

pet_dpeak <- lapply(pet_dpeak,dpeakFit,maxComp = g,nCore = mc)

pet_bs_loc <- file.path(bases[[2]],paste0("pet_bs_g",g,"_",1:2,".bed"))

for(i in 1:2){
  export(pet_dpeak[[i]],type = "bed",filename = pet_bs_loc[i])
}


## chip set pet binding sites

set_read_loc <- file.path(set_dir,read_files(set_dir,set,FALSE))

set_dpeak <- mapply(dpeakRead,set_peaks_loc,set_read_loc,
  MoreArgs = list(fileFormat = "bam",PET = FALSE,fragLen = fl,parallel = TRUE,nCore = mc),SIMPLIFY =FALSE)

set_dpeak <- lapply(set_dpeak,dpeakFit,maxComp = g,nCore = mc)

set_bs_loc <- file.path(bases[[3]],paste0("set_bs_g",g,"_",1:2,".bed"))

for(i in 1:2){
  export(set_dpeak[[i]],type = "bed",filename = set_bs_loc[i])
}

######################################################################################
