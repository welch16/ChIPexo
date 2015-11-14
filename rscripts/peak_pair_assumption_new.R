
rm(list = ls())

library(mosaics)
library(ChIPUtils)

library(dpeak)

exo_file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo/edsn931_042814_qc.sorted.bam"
pet_file1 <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET/edsn790_042814_qc.filter.bam"
set_file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET/edsn80_042814_qc.sorted.bam"

frag_len <- 150
bin_size <- 150

extra <- list()
extra[["exo"]] <- file.path("/p/keles/genome_data/EColi_U00096.2",
            c("mappability",
              "GC",
              "N"),"bin",
            c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
              "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
              "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))
extra[["pet"]] <- "/p/keles/ChIPexo/volume6/inputs/pet/edsn1369_042814_qc.filter.bam_bin150.txt"
extra[["set"]] <- "/p/keles/ChIPexo/volume6/inputs/set/edsn101_042814_qc.sorted.bam_fragL150_bin150.txt"


check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)


out_dir <- "/p/keles/ChIPexo/volume6/peak_pair"

### create bins


bin_flag <- FALSE

bin_dir <- file.path(out_dir,"bins")
check_create(bin_dir)


if(bin_flag){
  constructBins(exo_file,fileFormat = "bam",outfileLoc = bin_dir,PET = FALSE,
                fragLen = frag_len,binSize = bin_size)
  constructBins(pet_file1,fileFormat = "bam",outfileLoc = bin_dir, PET =  TRUE,
               binSize = bin_size)
  constructBins(pet_file2,fileFormat = "bam",outfileLoc = bin_dir, PET =  TRUE,
               binSize = bin_size)
  constructBins(set_file,fileFormat = "bam",outfileLoc = bin_dir, PET =  FALSE,
               fragLen = frag_len , binSize = bin_size)  
}

### call peaks



peak_flag <- FALSE

if(peak_flag){
  exo_bin <- file.path(bin_dir,paste0(basename(exo_file),"_fragL",frag_len,"_bin",bin_size,".txt"))
  pet_bin1 <- file.path(bin_dir,paste0(basename(pet_file1),"_bin",bin_size,".txt"))
  pet_bin2 <- file.path(bin_dir,paste0(basename(pet_file2),"_bin",bin_size,".txt"))
  set_bin <- file.path(bin_dir,paste0(basename(set_file),"_fragL",frag_len,"_bin",bin_size,".txt"))
  
  exo_bins <- readBins(c("chip","M","GC","N"),c(exo_bin,extra[["exo"]]))
  pet_bins1 <- readBins(c("chip","input"),c(pet_bin1,extra[["pet"]]))
  pet_bins2 <- readBins(c("chip","input"),c(pet_bin2,extra[["pet"]]))
  set_bins <- readBins(c("chip","input"),c(set_bin,extra[["set"]]))

  
  exo_fit <- mosaicsFit(exo_bins)
  pet_fit1 <- mosaicsFit(pet_bins1)
  pet_fit2 <- mosaicsFit(pet_bins2)
  set_fit <- mosaicsFit(set_bins)

  pdf(file = "figs/sig70/gof_plots_sig70.pdf")
  plot(exo_fit)
  plot(pet_fit1)
  plot(pet_fit2)
  plot(set_fit)
  dev.off()


  get_peaks <- function(fit,fdr,bin_size,thresh)
  {
    if(fit@bic2S <= fit@bic1S){
      peak <- mosaicsPeak(fit , signalModel = "2S", FDR = fdr, maxgap = bin_size ,thres = thresh)
    }else{
      peak <- mosaicsPeak(fit , signalModel = "1S", FDR = fdr, maxgap = bin_size ,thres = thresh)
    }
    show(peak)
    return(data.table(peak@peakList))
  }
  
  fdr <- .05
  thresh <- 10

  exo_peak <- get_peaks(exo_fit,fdr,bin_size,thresh)
  pet_peak1 <- get_peaks(pet_fit1,fdr,bin_size,thresh)
  pet_peak2 <- get_peaks(pet_fit2,fdr,bin_size,thresh)
  set_peak <- get_peaks(set_fit,fdr,bin_size,thresh)

  peak_dir <- file.path(out_dir,"peaks")
  check_create(peak_dir)
  exo_peak_file <- file.path(peak_dir,gsub("_qc.sorted.bam","_peaks.txt",basename(exo_file)))
  pet_peak_file1 <-  file.path(peak_dir,gsub("_qc.filter.bam", "_filter_peaks.txt",basename(pet_file1)))
  pet_peak_file2 <- file.path(peak_dir, gsub("_qc.sorted.bam", "_all_peaks.txt",basename(pet_file2)))
  set_peak_file <-  file.path(peak_dir,gsub("_qc.sorted.bam", "_peaks.txt",basename(set_file)))

  

  write.table(exo_peak,file = exo_peak_file, sep = "\t",row.names = FALSE,quote = FALSE,col.names=FALSE)
  write.table(pet_peak1, file  = pet_peak_file1, sep = "\t",row.names = FALSE,quote = FALSE,col.names=FALSE)
  write.table(pet_peak2, file  = pet_peak_file2, sep = "\t",row.names = FALSE,quote = FALSE,col.names=FALSE)
  write.table(set_peak , file = set_peak_file, sep = "\t",row.names = FALSE,quote = FALSE,col.names=FALSE)
}

  peak_dir <- file.path(out_dir,"peaks")

  ## exo_peak_file <- file.path(peak_dir,gsub("_qc.sorted.bam","_peaks.txt",basename(exo_file)))
  pet_peak_file1 <-  file.path(peak_dir,gsub("_qc.filter.bam", "_filter_peaks.txt",basename(pet_file1)))
  ## pet_peak_file2 <- file.path(peak_dir, gsub("_qc.sorted.bam", "_all_peaks.txt",basename(pet_file2)))
  ## set_peak_file <-  file.path(peak_dir,gsub("_qc.sorted.bam", "_peaks.txt",basename(set_file)))


## estimate binding events

bs_dir <- file.path(out_dir,"binding")
check_create(bs_dir)

bind_flag <- TRUE

## load_all("~/Desktop/Docs/Code/dpeak")



if(bind_flag){

  ## exo_dp <- dpeakRead(peakfile = exo_peak_file,readfile = exo_file , fileFormat = "bam",
  ##                 PET = FALSE, fragLen = frag_len,parallel = TRUE,nCore = 18)
  ## exo_fits <- list()
  ## exo_fits[["common"]] <- dpeakFit(exo_dp,nCore = 18) 
  ## exo_fits[["separate"]] <- dpeakFit(exo_dp,estDeltaSigma = "separate",nCore = 18)
  ## save( exo_fits, file = file.path(bs_dir, "dpeakFit_exo_sig70aerobic.RData"))

  pet_dp1 <- dpeakRead(peakfile = pet_peak_file1,readfile = pet_file1 , fileFormat = "bam",
                  PET = TRUE, fragLen = frag_len,parallel = TRUE,nCore = 24)   
  pet_fits1 <- list()
  pet_fits1[["common"]] <- dpeakFit(pet_dp1,nCore = 24) 
  ## pet_fits1[["separate"]] <- dpeakFit(pet_dp1,estDeltaSigma = "separate",nCore = 18)
  save( pet_fits1, file = file.path(bs_dir, "dpeakFit_pet3_sig70aerobic.RData"))


  ## pet_dp2 <- dpeakRead(peakfile = pet_peak_file2,readfile = pet_file2 , fileFormat = "bam",
  ##                 PET = TRUE, fragLen = frag_len,parallel = TRUE,nCore = 18)
  ## pet_fits2 <- list()
  ## pet_fits2[["common"]] <- dpeakFit(pet_dp2,nCore = 18) 
  ## pet_fits2[["separate"]] <- dpeakFit(pet_dp2,estDeltaSigma = "separate",nCore = 18)
  ## save( pet_fits2, file = file.path(bs_dir, "dpeakFit_pet2_sig70aerobic.RData"))
  
  ## set_fits <- list()
  ## set_dp <- dpeakRead(peakfile = set_peak_file,readfile = set_file , fileFormat = "bam",
  ##                 PET = FALSE, fragLen = frag_len,parallel = TRUE,nCore = 18)
  ## set_fits[["common"]] <- dpeakFit(set_dp,nCore = 18) 
  ## set_fits[["separate"]] <- dpeakFit(set_dp,estDeltaSigma = "separate",nCore = 18)    
  ## save( set_fits, file = file.path(bs_dir, "dpeakFit_set_sig70aerobic.RData"))

}

## load( file = file.path(bs_dir, "dpeakFit_sig70aerobic.RData"))

## sites1 <- data.table(read.table(file.path(bs_dir,"sig70_bs_common.bed"),skip = 1))
## sites2 <- data.table(read.table(file.path(bs_dir,"sig70_bs_separate.bed"),skip = 1))

## peaks[ , fsr := sapply(dp1@fragSet, function(x) mean(x[,3] == "F"))]


## library(ggplot2)
## library(scales)





  
