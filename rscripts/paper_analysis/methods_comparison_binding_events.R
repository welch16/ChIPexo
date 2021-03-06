
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(ChIPUtils)
library(gridExtra)
library(grid)
library(parallel)

######################################################################################

## Initial parameters

mc <- detectCores()
figs_dir <- "figs/for_paper"

dir.create(figs_dir,showWarnings = FALSE, recursive = TRUE)
base_dir <- "/p/keles/ChIPexo/volume6/K12"
FDR <- "FDR1"
edsn <- c("1311","1314","1317","1320","931","933")
max_dist <- 100
topM <- 200
mc <- detectCores()

### methods directories

peakzilla_dir <- file.path(base_dir,"other_methods/peakzilla")
mace_dir <- file.path(base_dir,"other_methods/mace")
gem_dir <- file.path(base_dir,"other_methods/gem")

######################################################################################

## get annotations -- gold standard

annot_dir <- file.path(base_dir,"annotations")
afiles <- list.files(annot_dir)
afiles <- afiles[grep("csv",afiles)]

annots <- read.csv(file.path(annot_dir,afiles),sep = "\t")
annots <- data.table(annots)
annots <- annots[grep("Sigma70",sigma.factor)]
annots <- annots[grep("TIM",description)]
annots <- annots[grep("RPP",description,invert = TRUE)]
annots <- annots[grep("AIPP",description,invert = TRUE)]
annots <- annots[grep("WHO",description,invert = TRUE)]
annots <- annots[grep("NTAS",description,invert = TRUE)]
annots <- annots[grep("ICA",description,invert = TRUE)]
annotDT <- annots

annots <- annots[,IRanges(start = coord,width = 1)]
annots <- sort(annots)

######################################################################################

## load method sites

### METHOD 1) load peakzilla sites

files <- list.files(peakzilla_dir)
files <- files[grep("tsv",files)]

peakzilla <- lapply(file.path(peakzilla_dir,files),read.table)
peakzilla <- lapply(peakzilla,data.table)
peakzilla <- lapply(peakzilla,function(x){
  setnames(x,names(x),c("Chromosome","Start","End","Name","Summit","Score","ChIP","Control",
                        "FoldEnrichment","DistributionScore","FDR"))
  x <- x[order(Start)]
  return(x)})
names(peakzilla) <- edsn

rm(files)

### METHOD 2) load mace sites

files <- list.files(mace_dir)
files <- files[grep("bed",files)]
files <- files[grep("border_pair.bed",files,fixed = TRUE)]

mace <- lapply(file.path(mace_dir,files),read.table)
mace <- lapply(mace, data.table)

mace <- lapply(mace,function(x){
  setnames(x,names(x),c("chrom","st","en","border","pval"))
  x[,Name := paste0(chrom,":",st,"-",en)]
  x[,site := mid(IRanges(start = st , end = en))]
  return(x)})
names(mace) <- edsn
rm(files)

### METHOD 3) load gem sites

files <- list.files(gem_dir,recursive = TRUE)
files <- files[grep(FDR,files)]
if(FDR == "FDR1"){
  files <- files[grep("FDR10",files,invert = TRUE)]
}
files <- files[grep("bed",files)]
files <- files[grep("insig",files,invert = TRUE)]
files <- files[grep("event",files)]

gem <- lapply(file.path(gem_dir,files),read.table,skip = 1)
gem <- lapply(gem,data.table)

gem <- lapply(gem,function(x){
  setnames(x,names(x),c("chr","st","en","name","signal"))
  x[,site := mid(IRanges(start = st , end = en))]
  return(x)})
names(gem) <- edsn

rm(files)

######################################################################################

## Get mosaics peaks and dpeak sites in two separate lists

### A) MOSAiCS peaks

files <- list.files(file.path(base_dir,"downstream"),recursive = TRUE)
files <- files[grep("peak",files)]
files <- files[grep("exo",files)]
files <- files[grep(FDR,files)]
if(FDR == "FDR1"){
  files <- files[grep("FDR10",files,invert =TRUE)]
}

mosaics <- lapply(file.path(base_dir,"downstream",files),read.table)
mosaics <- lapply(mosaics,data.table)
mosaics <- lapply(mosaics,function(x){
  setnames(x,names(x),c("chrID","peakStart","peakStop","peakSize","logAveP","logMinP",
                        "aveLogP","aveChipCount","maxChipCount","map","GC"))
  x})
names(mosaics) <- edsn
rm(files)

mosaics <- lapply(mosaics,function(x,topM)x[order(-aveChipCount)][1:topM],topM)
mosaics_ranges <- lapply(mosaics,function(x) x[,IRanges(start = peakStart, end = peakStop)])
names(mosaics_ranges) <- edsn

### B) dPeak sites
files <- list.files(file.path(base_dir,"downstream"),recursive = TRUE)
files <- files[grep("exo",files)]
files <- files[grep("sites",files)]
files <- files[grep(FDR,files)]
if(FDR == "FDR1"){
  files <- files[grep("FDR10",files,invert = TRUE)]
}
files_mot <- files[grep("motif",files)]
files <- files[grep("motif",files,invert = TRUE)]

dpeak <- lapply(file.path(base_dir,"downstream",files),read.table,header = TRUE)
dpeak <- lapply(dpeak,data.table)
names(dpeak) <- edsn
rm(files)

dpeak_mot <- lapply(file.path(base_dir,"downstream",files_mot),read.table,header = TRUE)
dpeak_mot <- lapply(dpeak_mot,data.table)
names(dpeak_mot) <- edsn
rm(files_mot)

dpeak_mot <- lapply(dpeak_mot,function(x){
  setnames(x,names(x),c("chrID","siteStart","siteEnd","siteName","siteStrength"))
  return(x)})


### C)dPeak sites with motif init

tab1 <- do.call(rbind,list(peakzilla = sapply(peakzilla,nrow),
                           mace = sapply(mace,nrow),
                           gem = sapply(gem,nrow),
                           mosaics = sapply(mosaics,nrow),
                           dpeak = sapply(dpeak,nrow),
                           dpeak_mot = sapply(dpeak_mot,nrow)))
tab1

######################################################################################

covert2IRanges <- function(DT,method)
{
  if(method == "peakzilla"){
    out <- DT[,IRanges(start = mid(IRanges(start = Start, end = End)),width = 1)]
  }else if(method == "mace"){
    out <- DT[,IRanges(start = mid(IRanges(start = st,end = en)),width = 1)]
  }else if(method == "gem"){
    out <- DT[,IRanges(start = mid(IRanges(start = st,end = en)),width = 1)]
  }else if(method == "dpeak"){
    out <- DT[,IRanges(start = mid(IRanges(start = siteStart,end = siteEnd)),width = 1)]
  }else{
    message(method , " is unknown, returning empty IRanges")
    out <- IRanges()
  }
  return(out)
}

## Convert files to sites

peakzilla_ranges <- lapply(peakzilla,covert2IRanges,"peakzilla")
mace_ranges <- lapply(mace,covert2IRanges,"mace")
gem_ranges <- lapply(gem,covert2IRanges,"gem")
dpeak_ranges <- lapply(dpeak,covert2IRanges,"dpeak")
dpeak_mot_ranges <- lapply(dpeak_mot,covert2IRanges,"dpeak")


######################################################################################

## Compare ranges
resolution <- function(annot,predictions,max_dist)
{
  if(length(predictions) == 0){
    out <- NA
  }else{
    out <- sapply(mid(annot),function(x,pred)min(abs(x - pred)),mid(predictions))
    if(length(predictions) < length(annot)){
      out <- sort(out)[1:length(predictions)]
    }
    out[out > max_dist] <- NA
  }
  return(out)
}

separate_sites <- function(peaks,sites,idx = NULL)
{
  out <- mclapply(peaks,function(x)subsetByOverlaps(sites,x),mc.cores = mc)
  if(!is.null(idx))out <- out[idx]
  return(out)
}

clean_list <- function(predictions,nms)
{
  nn1 <- nms[row %in% names(predictions)]
  names(predictions) <- nn1[,(peaks)]
  out <- mapply(function(x,y)data.table(reso = x,peak = y),
      predictions,names(predictions),
      SIMPLIFY = FALSE)
  out <- do.call(rbind,out)
  out <- out[!is.na(reso)]
  return(out)
}

compare_predictions <- function(peakzilla,mace,gem,dpeak,dpeak_mot,peaks,annots,max_dist)
{
  peaks <- split(peaks,1:length(peaks))
  annots_in_peaks <- separate_sites(peaks,annots)  
  idx <- which(sapply(annots_in_peaks,length) > 0)
  nms <- data.table(peaks =
   sapply(peaks,function(x)paste(start(x),end(x),sep = "-")),
   row = as.character(1:length(peaks)))
  annots_in_peaks <- annots_in_peaks[idx]
  peakzilla_in_peaks <- separate_sites(peaks,peakzilla,idx)
  mace_in_peaks <- separate_sites(peaks,mace,idx)
  gem_in_peaks <- separate_sites(peaks,gem,idx)
  dpeak_in_peaks <- separate_sites(peaks,dpeak,idx)
  dpeak_mot_in_peaks <- separate_sites(peaks,dpeak_mot,idx)
  
  peakzilla_reso <- mcmapply(resolution,annots_in_peaks,peakzilla_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  mace_reso <- mcmapply(resolution,annots_in_peaks,mace_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  gem_reso <- mcmapply(resolution,annots_in_peaks,gem_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  dpeak_reso <- mcmapply(resolution,annots_in_peaks,dpeak_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  dpeak_mot_reso <- mcmapply(resolution,annots_in_peaks,dpeak_mot_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)

  peakzilla_reso <- clean_list(peakzilla_reso,nms)[,method :="Peakzilla"]
  mace_reso <- clean_list(mace_reso,nms)[,method := "Mace"]
  gem_reso <- clean_list(gem_reso,nms)[,method :=  "Gem"]
  dpeak_reso <- clean_list(dpeak_reso,nms)[,method := "dPeak"]
  dpeak_mot_reso <- clean_list(dpeak_mot_reso,nms)[,method := "dPeak (motif)"]
  out <- do.call(rbind,list(peakzilla_reso,mace_reso,gem_reso,dpeak_reso,dpeak_mot_reso))
  return(out)
}

reso <- mapply(compare_predictions,
               peakzilla_ranges,
               mace_ranges,
               gem_ranges,
               dpeak_ranges,
               dpeak_mot_ranges,
               mosaics_ranges,
               MoreArgs = list(annots,max_dist),SIMPLIFY = FALSE)
reso <- mapply(function(x,y)copy(x)[,dataset := y],reso,edsn,SIMPLIFY = FALSE)
reso <- do.call(rbind,reso)

save(reso,file = "data/Resolution_method_with_peaks.RData")

pdf(file = "figs/gem_dpeak_reso/nr_sites.pdf",height = 12,width = 4)
ggplot(melt(dcast.data.table(reso[dataset == 931], value.var = "reso",
     formula = peak + dataset ~ method,
     fun.aggregate = length)
),aes( variable,peak , fill = factor(value)))+geom_tile()+
  theme_bw()+
  theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
  scale_fill_brewer(palette = "Pastel1",name = "nsites")
dev.off()

reso2 <- dcast.data.table(reso[dataset == 931], value.var = "reso",
     formula = peak + dataset ~ method,
     fun.aggregate = length)
aux <- strsplit(reso2[,(peak)],"-")
dt <- data.table(sq = "U00096",
   start = sapply(aux,function(x)as.numeric(x[1])),
   end = sapply(aux,function(x)as.numeric(x[2])))

dt[start == 0, start := 1]

write.table(dt,
 file = "/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test/peak_to_check_dpeak.txt",
 quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

u <- dt[, paste0(sq,":",start,"-",end)]
write.table(u,
 file = "/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test/peak_to_check_gem.txt",
 quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
