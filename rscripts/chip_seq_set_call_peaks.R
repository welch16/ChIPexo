
rm(list = ls())

library(mosaics)
library(parallel)
library(data.table)

mc <- 12


bin_dir <- "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_SET/bins"
input_dir <- "/p/keles/genome_data/EColi_U00096.2/"
figs_dir <- "figs/mosaics"

binSize <- 150
fragLen <- 150

files <- list.files(bin_dir) 

## GC_file <- file.path(input_dir,
##   "GC/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_GC_fragL150_bin150.txt")
## map_file <- file.path(input_dir,
##   "mappability/bin/fragL150_bin150_32mer_single/E.coli_K-12_MG1655_GENOME_fragL150_bin150.txt")
## N_file <- file.path(input_dir,
##   "N/bin/fragL150_bin150/E.coli_K-12_MG1655_GENOME_N_fragL150_bin150.txt")

opt <- c("chip","input")

read_bins_wrap <- function(file,input_file,opt)
{
  readBins(opt,c(file, input_file))
}


## files[10] is a beta sample, which have some issues, (that's way they remapped them to beta'_f
## files[13] is the input for chip_exo...


bins <- mclapply(file.path(bin_dir,files[-1]),read_bins_wrap, file.path(bin_dir,files[1]),opt,mc.cores = mc)

## exploratory plots

pdf(file = file.path(figs_dir,"chip_seq_set_bins.pdf"))
u <- lapply(bins,plot)
dev.off()


fits <- mclapply(bins,mosaicsFit,mc.cores = mc)

## from observing the the gof plots, we can see that only one that is kinda bed is the one for the
## edsn932 file, which corresponds to sigma_s - exponential growth - replicate 1

pdf(file = file.path(figs_dir,"chip_seq_set_gof.pdf"))
u <- lapply(fits,plot)
dev.off()
    
FDR <- .25
thresh <- 10

mosaics_peak_wrap <- function(fit,FDR,binsize,thres)
{
  if ( fit@bic2S <= fit@bic1S ) {
    peak <- mosaicsPeak( fit, signalModel="2S", FDR=FDR, maxgap=binsize, thres=thres )
  } else {
    peak <- mosaicsPeak( fit, signalModel="1S", FDR=FDR, maxgap=binsize, thres=thres )
  }
  return(peak)
}

peaks <- mclapply(fits,mosaics_peak_wrap,FDR,binSize,thresh,mc.cores = mc)

dt_list <- lapply(peaks,function(x)data.table(x@peakList))

out_dir <- "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_SET/peaks"

out_files <- gsub(".filter.bam_fragL150_bin150.txt","_peaks.txt",files)[-1]


for(i in 1:length(out_files)){
  write.table(dt_list[[i]],file = file.path(out_dir,out_files[i]),col.names =  FALSE,row.names = FALSE,quote = FALSE)
}
