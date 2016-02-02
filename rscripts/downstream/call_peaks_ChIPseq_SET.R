rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(parallel)

frag_len <- 150
bin_size <- 150
fdr <- .05

out_dir <- "/p/keles/ChIPexo/volume6/resolution/ChIPseq_SET"
in_dir <- "/p/keles/ChIPexo/volume7/Landick/ChIPseq_SET/rif_treatment"

files <- c("edsn1396_Sig70.sort.bam","edsn1398_Sig70.sort.bam",
           "edsn1400_Sig70.sort.bam","edsn1402_Sig70.sort.bam")

binfiles <- paste0(files,"_fragL",frag_len,"_bin",bin_size,".txt")

bin_dir <- file.path(out_dir,"bins")

extra <- "/p/keles/ChIPexo/volume6/resolution/inputs/ChIPseq_SET/edsn1369_Input.sort.bam_fragL150_bin150.txt"

bins <- mclapply(file.path(bin_dir,binfiles),
  function(x)readBins(type = c("chip","input"),
                      fileName = c(x,extra)),mc.cores = 4)

fits <- mclapply(bins,mosaicsFit,mc.cores = 4)

## lapply(fits,plot)
## dev.off()

FF <- paste0("FDR",fdr*100)
out_dir <- file.path(out_dir,"peaks",FF)
dir.create(out_dir,recursive = TRUE , showWarnings = FALSE)

peakfiles <- gsub(".sort.bam","_peaks.txt",files)

peaks <- mclapply(fits,function(x){
  model <- ifelse(x@bic1S < x@bic2S,"1S","2S")
  mp <- mosaicsPeak(x,signalModel = model,thres = 300)
  out <- data.table(print(mp))
  return(out)},mc.cores = 4)

mapply(write.table,peaks,file.path(out_dir,peakfiles),
       MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE,
         sep = "\t"))
