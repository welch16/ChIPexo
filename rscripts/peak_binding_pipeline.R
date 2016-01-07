rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(dpeak)

file <- "/p/keles/ChIPexo/volume7/rif_treatment/edsn1320_Sig70.sort.bam"
out_dir <- "/p/keles/ChIPexo/volume6/bowtie_alignments/example"

frag_len <- 150
bin_size <- 150
fdr <- .05
thresh <- 10
mc <- 24
pet <- FALSE

if(pet){

}else{
  extra <- file.path("/p/keles/genome_data/EColi_U00096.2",
            c("mappability",
              "GC",
              "N"),"bin",
            c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
              "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
              "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))
}

  
binfile <- paste0(basename(file),"_fragL",frag_len,"_bin",bin_size,".txt")
peakfile <- gsub(".sort.bam","_peaks.txt",basename(file))
bindingfile <- gsub(".sort.bam", "_binding.RData",basename(file))


## extra[["pet"]] <- "/p/keles/ChIPexo/volume6/inputs/pet/edsn1369_042814_qc.filter.bam_bin150.txt"
## extra[["set"]] <- "/p/keles/ChIPexo/volume6/inputs/set/edsn101_042814_qc.sorted.bam_fragL150_bin150.txt"

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

##### flags

bin_flag <- TRUE
peak_flag <- TRUE
binding_ev_flag <- TRUE

check_create(out_dir)

bin_dir <- file.path(out_dir,"bins")
check_create(bin_dir)

### bins

if(bin_flag){
  constructBins(file,fileFormat = "bam",outfileLoc = bin_dir,PET = pet,
                fragLen = frag_len,binSize = bin_size)
}

### call peaks

get_peaks <- function(fit,fdr,bin_size,thresh)
{
  if(fit@bic2S <= fit@bic1S){
    peak <- mosaicsPeak(fit , signalModel = "2S", FDR = fdr, maxgap = bin_size ,thres = thresh)
  }else{
    peak <- mosaicsPeak(fit , signalModel = "1S", FDR = fdr, maxgap = bin_size ,thres = thresh)
  }
  show(peak)
  return(peak)
}


save_fix <- function(peak,file)
{
  dt <- data.table(peak@peakList)
  write.table(dt , file = file , quote = FALSE,col.names = FALSE,row.names = FALSE)
}




if(peak_flag){
  binfile <- file.path(bin_dir,binfile)
  
  bins <- readBins(c("chip","M","GC","N"),c(binfile,extra))

  fit <- mosaicsFit(bins)
  plot(fit)
  ## wanna make a ggplot version

  peaks <- get_peaks(fit,fdr,bin_size,thresh)
  
  peak_dir <- file.path(out_dir,"peaks")
  check_create(peak_dir)
  peakfile <- file.path(peak_dir,peakfile)
  save_fix(peaks,peakfile)

}

bs_dir <- file.path(out_dir,"binding")
check_create(bs_dir)

if(binding_ev_flag){  
  dp <- dpeakRead(peakfile = peakfile,readfile = file , fileFormat = "bam",
                  PET = FALSE, fragLen = frag_len,parallel = TRUE,nCore = mc)
  fit <- dpeakFit(dp,nCore = mc)
  save(fit, file = file.path(bs_dir,bindingfile))

}

