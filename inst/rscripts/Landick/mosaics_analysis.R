
# only sig70

rm(list = ls())

library(mosaics)
library(data.table)
library(ggplot2)
library(parallel)
library(RColorBrewer)
library(hexbin)
library(scales)

filesdir = "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET/sam"
outdir = "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET/bins"
figsdir = "/p/keles/ChIPexo/volume3/ChIPexo/inst/figs/Landick/mosaics"

# sig 70
## edsn930_042814_qc.sorted.bed input
## edsn931_042814_qc.sorted.bed exp-rep1
## edsn933_042814_qc.sorted.bed exp-rep2
## edsn935_042814_qc.sorted.bed stat-rep1
## edsn937_042814_qc.sorted.bed stat-rep2
## edsn1311_042814_qc.sorted.bed rif0 rep1
## edsn1314_042814_qc.sorted.bed rif20 rep1
## edsn1317_042814_qc.sorted.bed rif0 rep2
## edsn1320_042814_qc.sorted.bed rif20 rep2

# sig 70
## edsn1369_042814_qc.sorted.bed seq-input
## edsn1396_042814_qc.sorted.bed rif0 rep1
## edsn1398_042814_qc.sorted.bed rif20 rep1
## edsn1400_042814_qc.sorted.bed rif0 rep2
## edsn1402_042814_qc.sorted.bed rif20 rep2



# Construct bins
files = c(
  "edsn1369_042814_qc.sorted.sam",
  "edsn1396_042814_qc.sorted.sam",
  "edsn1398_042814_qc.sorted.sam",
  "edsn1400_042814_qc.sorted.sam",
  "edsn1402_042814_qc.sorted.sam")

infiles = file.path(filesdir,files)

fragLen = 150
binSize = 150

## fragLen = 100 or 150 
## binSize = 100 or 150


# the files are named as the same regardless of capping
doBins = TRUE
if(doBins)
{
  u =mclapply(infiles, function(x){
    constructBins(infile = x,
      fileFormat = "sam",
      outfileLoc = outdir,
      byChr = FALSE,
      useChrfile = TRUE,
      chrfile = file.path(filesdir,"ecoli_size.txt"),
      PET = TRUE,
      fragLen = fragLen,
      binSize = binSize,
      capping = 0,
      perl = "perl")},mc.cores = 5)                
}

# enrichment plots

# read data
binfiles = file.path(outdir,
  paste0(files,"_bin",binSize,".txt"))

bins = mclapply(binfiles,function(x){
  dt = data.table(read.table(x,stringsAsFactors=FALSE))
  setnames(dt,names(dt),c("chr","coord","tagCount"))
  return(dt)},mc.cores = 9)
names(bins) =files

enrichment_plot <- function(bin_input,bin_chip,nrbins)
{
  rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r = rf(16)
  # enrichment plots, input vs chip hexbin plot
  dt = data.table(input = bin_input[,.(tagCount)][[1]],
    chip = bin_chip[,.(tagCount)][[1]])
  p = ggplot(dt,aes(input,chip))+stat_binhex(bins = nrbins)+
    scale_fill_gradientn(colours =r,trans='log10')+
    theme(legend.position = "top")
  return(p)
}


lapply(bins[-1],function(x){
  x11()
  print(enrichment_plot(bins[[1]],x,70))
})

# in all samples there is a vertical arm close to the chip axis, and therefore there is
# some enrichment in all samples

# this plots show two strong arms of bins close to both the input and
# chip axes. this suggests that as expected there is some regions with
# pcr artifacts. (this follows from the fact that this hexbin plots look very similar to
# the ones from the University of Nebraska's data sets)

# possibly this data sets are gonna work fine to call peaks


mos_fit <- function(chip,input,d)
{
  bin_calc = readBins(type = c("chip","input"),fileName = c(chip,input))
  fit = mosaicsFit(bin_calc,analysisType = "IO",bgEst = "automatic",d=d)
  return(fit)
}

d = .75
mosaics_fit_list = mclapply(binfiles[-1],function(x)mos_fit(x,binfiles[1],d),mc.cores = 9)
# length = 8 , 4 of the stat-vs-exp experiment and 4 of the rif-treatment experiment

## to plot gof
lapply(mosaics_fit_list,function(x){x11();plot(x)})


## parameters
## fdr = c(0.05,0.01,0.001,0.0001)
## maxgap = seq(100,1e3,by=1e2)
## threshold = seq(3,30,by = 3)


# we consider using a very conservative fdr = .1m low theshold to include low quality regions
# and we are merging by the binSize to only merge adjacent bins.
thresh = 1
FDR =.1
maxgap = binSize

mosaics_peaks_list = mclapply(mosaics_fit_list,function(x){
  mosaicsPeak(x,signalModel = "2S",thres = thresh,
    FDR = FDR,maxgap = maxgap)
},mc.cores =8)
names(mosaics_peaks_list) = files[-1]


peaks_list = lapply(mosaics_peaks_list,function(x)data.table(x@peakList))

save(peaks_list,file = file.path("/p/keles/ChIPexo/volume3/Analysis/Landick","peaks.RData"))

    ## Created the script mosaics_analysis to call a conservative collection
    ## of peaks for all Landick's data sets.
    
    ## All peaks were called using the "Input-Only" fit, using d = 0.75. This
    ## parameter was selected by fiting all data sets with a mosaics model,
    ## and selecting the best BIC from d in (.1,.5,.75,1)
    
    ## Then we selected the optimal d as the one that was the best for most datasets.
    
    ## The peaks were called using
    ## thres = 1
    ## FDR = .1
    ## maxgap = binSize
    
    ## to obtain a conservative list of peaks


