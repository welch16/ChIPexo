rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(parallel)

frag_len <- 150
bin_size <- 150
fdr <- .01


out_dir <- "/p/keles/ChIPexo/volume6/resolution/ChIPexo"
in_dir <- "/p/keles/ChIPexo/volume7/Landick/ChIPexo/rif_treatment"

files <- c("edsn1311_Sig70.sort.bam","edsn1314_Sig70.sort.bam",
           "edsn1317_Sig70.sort.bam","edsn1320_Sig70.sort.bam")

binfiles <- paste0(files,"_fragL",frag_len,"_bin",bin_size,".txt")

bin_dir <- file.path(out_dir,"bins")

extra <- file.path("/p/keles/genome_data/EColi_U00096.2",
          c("mappability","GC","N"),
          "bin",
          c("fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
            "fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
            "fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))

bins <- mclapply(file.path(bin_dir,binfiles),
  function(x)readBins(type = c("chip","M","GC","N"),
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
  mp <- mosaicsPeak(x,signalModel = model,thres = 100)
  out <- data.table(print(mp))
  return(out)},mc.cores = 4)

mapply(write.table,peaks,file.path(out_dir,peakfiles),
       MoreArgs = list(quote = FALSE,row.names = FALSE,sep = "\t"))
