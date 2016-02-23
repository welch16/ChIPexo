rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(parallel)

frag_len <- 150
bin_size <- 150

in_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo"

files <- c("edsn1311_Sig70.sort.bam","edsn1314_Sig70.sort.bam",
           "edsn1317_Sig70.sort.bam","edsn1320_Sig70.sort.bam",
           "edsn931_Sig70.sort.bam","edsn933_Sig70.sort.bam")

binfiles <- paste0(files,"_fragL",frag_len,"_bin",bin_size,".txt")

bin_dir <- file.path(in_dir,"bins")

extra <- file.path("/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/peak_inputs",
          c("mappability","GC","N"),
          c("E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt",
            "E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt",
            "E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt"))

bins <- mclapply(file.path(bin_dir,binfiles),
  function(x)readBins(type = c("chip","M","GC","N"),
                      fileName = c(x,extra)),mc.cores = 6)

fits <- mclapply(bins,mosaicsFit,mc.cores = 6)

pdf(file = "figs/Sig70_mosaics_GOF.pdf")
lapply(fits,plot)
dev.off()

fdr <- .1
peaks <- mclapply(fits,function(x){
  model <- ifelse(x@bic1S < x@bic2S,"1S","2S")
  mp <- mosaicsPeak(x,signalModel = model,thres = 100,maxgap = 300)
  out <- data.table(print(mp))
  return(out)},mc.cores = 6)

FF <- paste0("FDR",fdr*100)
out_dir <- file.path(in_dir,"peaks",FF)
dir.create(out_dir,recursive = TRUE , showWarnings = FALSE)

peakfiles <- gsub(".sort.bam","_peaks.txt",files)

mapply(write.table,peaks,file.path(out_dir,peakfiles),
       MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t"))
