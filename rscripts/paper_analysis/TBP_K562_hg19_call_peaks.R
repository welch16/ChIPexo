
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
## library(dpeak)
library(data.table)
library(parallel)

indir = "/p/keles/ChIPexo/volume4"
files = list.files(indir,full.names = TRUE,recursive = TRUE)
files = files[grep("TBP",files)]
files = files[grep("bam",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("subsam",files,invert = TRUE)]
files = files[grep("chipseq",files,invert = TRUE)]


outdir = "/p/keles/ChIPexo/volume4/venters_data/my_peaks"

### construct bins
fragLen <- 200
binSize <- 200

options(mc.cores = 22)


## A = lapply(files,
##   constructBins,fileFormat = "bam",
##               outfileLoc = file.path(outdir,"bins"),
##               fragLen = fragLen ,binSize = binSize)

binfiles = list.files(outdir,full.names = TRUE,recursive = TRUE)
binfiles = binfiles[c(grep("bam",binfiles),
                      grep("all",binfiles))]

binfiles1 = binfiles[seq_len(5)]
binfiles2 = binfiles[-seq_len(5)]
names(binfiles2) = c("GC","M","N")

bins = lapply(binfiles1,
                function(x)readBins(type = c("chip","M","GC","N"),
                                    fileName = c(x,binfiles2[["M"]],
                                                 binfiles2[["GC"]],
                                                 binfiles2[["N"]]),
                                    parallel = TRUE,nCore = getOption("mc.cores")))

load(file = file.path(outdir,"TBP_fits.RData"))

fits = list()
fits[[1]] = mosaicsFit(bins[[1]],parallel = TRUE,nCore = getOption("mc.cores"))
fits[[2]] = mosaicsFit(bins[[2]],parallel = TRUE,nCore = getOption("mc.cores"))
## fits[[3]] = mosaicsFit(bins[[3]],parallel = TRUE,nCore = getOption("mc.cores"))
## fits[[4]] = mosaicsFit(bins[[4]],parallel = TRUE,nCore = getOption("mc.cores"))
## fits[[5]] = mosaicsFit(bins[[5]],parallel = TRUE,nCore = getOption("mc.cores"))



save(fits,file = file.path(outdir,"TBP_fits.RData"))

## pdf(file = file.path(outdir,"GOF_plots.pdf"))
## plot(fits[[1]])
## lapply(fits,plot)
## dev.off()

load(file.path(outdir,"TBP_fits.RData"))

call_peaks <- function(fit,fdr,thres,mg)
{
   
  mp <- mosaicsPeak(fit,signalModel = "2S",thres = thres,maxgap = mg)
  peaks <- data.table(print(mp))
  
  return(peaks)
}

fdr <- 0.05
thres <- 100
mg <- 200

peaks <- lapply(fits,call_peaks,fdr,thres,mg)

library(magrittr)

nms = basename(binfiles1) %>% strsplit("\\.") %>%
    sapply(function(x)x[1])

mapply(write.table,peaks,
       file = file.path(outdir,paste0(nms[1:2], "_peaks.txt")),
       MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

# chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC 





