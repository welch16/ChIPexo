
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
FDR <- "FDR5"
edsn <- c("1311","1314","1317","1320","931","933")
max_dist <- 40
topM <- 500
strength <- 4e3
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

### B) dPeak sites
files <- list.files(file.path(base_dir,"downstream"),recursive = TRUE)
files <- files[grep("exo",files)]
files <- files[grep("sites",files)]
files <- files[grep(FDR,files)]
files <- files[grep("motif",files,invert = TRUE)]

dpeak <- lapply(file.path(base_dir,"downstream",files),read.table,header = TRUE)
dpeak <- lapply(dpeak,data.table)
dpeak <- lapply(dpeak,function(x)x[siteStrength > strength])
names(dpeak) <- edsn
rm(files)

tab1 <- do.call(rbind,list(peakzilla = sapply(peakzilla,nrow),
                           mace = sapply(mace,nrow),
                           gem = sapply(gem,nrow),
                           mosaics = sapply(mosaics,nrow),
                           dpeak = sapply(dpeak,nrow)))
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

clean_list <- function(predictions)
{
  out <- do.call(c,predictions)
  out <- out[!is.na(out)]
  return(out)
}

compare_predictions <- function(peakzilla,mace,gem,dpeak,peaks,annots,max_dist)
{
  peaks <- split(peaks,1:length(peaks))
  annots_in_peaks <- separate_sites(peaks,annots)
  idx <- which(sapply(annots_in_peaks,length) > 0)
  annots_in_peaks <- annots_in_peaks[idx]
  peakzilla_in_peaks <- separate_sites(peaks,peakzilla,idx)
  mace_in_peaks <- separate_sites(peaks,mace,idx)
  gem_in_peaks <- separate_sites(peaks,gem,idx)
  dpeak_in_peaks <- separate_sites(peaks,dpeak,idx)
  
  peakzilla_reso <- mcmapply(resolution,annots_in_peaks,peakzilla_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  mace_reso <- mcmapply(resolution,annots_in_peaks,mace_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  gem_reso <- mcmapply(resolution,annots_in_peaks,gem_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  dpeak_reso <- mcmapply(resolution,annots_in_peaks,dpeak_in_peaks,
    MoreArgs = list(max_dist),mc.cores = mc,SIMPLIFY = FALSE)
  peakzilla_reso <- data.table(method = "Peakzilla",reso = clean_list(peakzilla_reso))
  mace_reso <- data.table(method = "Mace",reso = clean_list(mace_reso))
  gem_reso <- data.table(method = "Gem",reso = clean_list(gem_reso))
  dpeak_reso <- data.table(method = "dPeak",reso = clean_list(dpeak_reso))

  out <- do.call(rbind,list(peakzilla_reso,mace_reso,gem_reso,dpeak_reso))
  return(out)
}

reso <- mapply(compare_predictions,
               peakzilla_ranges,
               mace_ranges,
               gem_ranges,
               dpeak_ranges,
               mosaics_ranges,
               MoreArgs = list(annots,max_dist),SIMPLIFY = FALSE)
reso <- mapply(function(x,y)copy(x)[,dataset := y],reso,edsn,SIMPLIFY = FALSE)
reso <- do.call(rbind,reso)

## pdf(file = file.path("figs/for_paper",paste0("sig70_methods_comparison_",FDR,"_topM",topM,".pdf")))
## ggplot(reso,aes_string("method","reso",fill = "method"))+geom_boxplot(outlier.shape = NA)+
##   scale_fill_brewer(palette = "Pastel1")+facet_wrap(~ dataset,nrow = 2 )+
##   theme_bw()+theme(legend.position = "none",axis.text.x = element_text(angle = 30))+
##   coord_cartesian(ylim = c(-max_dist*.05,max_dist*1.05))+
##   xlab("Methods")+ylab("Resolution")
## dev.off()

## pdf(file = file.path("figs/for_paper",paste0("sig70_methods_comparison_",FDR,"_topM",topM,".pdf")))

pdf(file = file.path("figs/for_paper",paste0("methods_comp_",FDR,"_topM",topM,".pdf")),height = 4,width = 4)
ggplot(reso[dataset %in% c("931")],aes_string("method","reso",fill = "method"))+geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 30),
        plot.title = element_text(hjust = 0))+
  xlab("Methods")+ylab("Resolution")+ylim(0,30)+ggtitle("A")
ggplot(reso[dataset %in% c("933")],aes_string("method","reso",fill = "method"))+geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 30),
        plot.title = element_text(hjust = 0))+
  xlab("Methods")+ylab("Resolution")+ylim(0,30)+ggtitle("B")
dev.off()


  ## geom_jitter(aes(colour = method),size = 1.2)+scale_color_brewer(palette = "Set1")+
tab1


library(dplyr)


my_reso = reso %>% filter(dataset %in% c("931","933")) %>%
  mutate(replicate = ifelse(dataset == "931","Rep 1","Rep 2")) %>%
  group_by(method)


pdf(file = "figs/for_paper/methods_comp_repl_together.pdf",height = 4, width = 8)
ggplot(my_reso,aes(method,reso,fill = method))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid( . ~ replicate)+
  scale_fill_brewer(palette = "Pastel2")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylim(0,30)
dev.off() 
  


