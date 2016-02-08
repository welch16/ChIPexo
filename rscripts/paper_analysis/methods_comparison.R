
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(ChIPUtils)
library(gridExtra)
library(grid)

######################################################################################

## Initial parameters

mc <- detectCores()
figs_dir <- "figs/condition/rif_treatment_reso"

dir.create(figs_dir,showWarnings = FALSE, recursive = TRUE)

base_dir <- "/p/keles/ChIPexo/volume6/resolution"


######################################################################################

## get promoter locations table

prom_dir <- "/p/keles/ChIPexo/volume6/gold_standard/Landick"
pfiles <- list.files(prom_dir)
pfiles <- pfiles[grep(".bed",pfiles)]
promoters <- lapply(file.path(prom_dir,pfiles),read.table)
promoters <- lapply(promoters,data.table)

promoters[[1]][,str := "R"]
promoters[[2]][,str := "F"]

sites <- do.call(rbind,promoters)
sites[, id := 1:nrow(sites)]
setnames(sites,names(sites),c("seqnames","start","end","strand","id"))

sites <- sites[,ifelse(strand == "R",end,start)]
sites <- IRanges(start = sites,width = 1)

######################################################################################

## get mosaics peaks
peakfiles <- list.files(base_dir,recursive = TRUE)
peakfiles <- peakfiles[grep("peak",peakfiles)]
peakfiles <- peakfiles[grep("exo",peakfiles)]
peakfiles <- peakfiles[grep("FDR5",peakfiles)]

mosaics_peaks <- lapply(file.path(base_dir,peakfiles),function(x)data.table(read.table(x)[,1:3]))
mosaics_peaks <- lapply(mosaics_peaks,
  function(x){
    setnames(x,names(x),c("seqnames","start","end"))
    x <- ranges(dt2gr(x))
    return(x)})
names(mosaics_peaks) <- basename(peakfiles)

######################################################################################

## sites to consider

sites_mosaics <- lapply(mosaics_peaks,function(x,sites)
  sort(reduce(subsetByOverlaps(sites,x))),sites)

######################################################################################


## peakzilla analysis

dr <- "inst/peakzilla_analysis"
files <- list.files(dr)
files <- files[grep("tsv",files)]

peakzilla <- lapply(file.path(dr,files),read.table)
peakzilla <- lapply(peakzilla,data.table)
peakzilla <- lapply(peakzilla,function(x){
  setnames(x,names(x),c("Chromosome","Start","End","Name","Summit","Score","ChIP","Control","FoldEnrichment","DistributionScore","FDR"))
  x <- x[order(Start)]
  return(x)})
sapply(peakzilla,nrow)

peakzilla <- mapply(function(mosaics,peakzilla){
  idx <- countOverlaps(IRanges(start = peakzilla[,(Start)],end = peakzilla[,(End)]),
                       mosaics) > 0
  out <- peakzilla[idx]
  return(out)},mosaics_peaks,peakzilla,SIMPLIFY = FALSE)
sapply(peakzilla,nrow)

peakzilla <- mapply(function(peakzilla,sites){
  pranges <- IRanges(start = peakzilla[,(Start)],end = peakzilla[,(End)])
#  ov <- findOverlaps(pranges,sites)  
  reso <- peakzilla[,.(reso = min(abs(Summit - start(sites)))),by = Name]
  setkey(reso,Name)
#  reso <- reso[as.character(peakzilla[queryHits(ov),(Name)])]
  return(reso)
},peakzilla,sites_mosaics,SIMPLIFY = FALSE)

peakzilla <- mapply(function(x,y)copy(x)[ , edsn := y],peakzilla,files,SIMPLIFY = FALSE)
peakzilla <- do.call(rbind,peakzilla)
peakzilla[,edsn := gsub("_Sig70_peaks.tsv","",edsn)]

peakzilla[, repl := 0]
peakzilla[, rif := ""]


peakzilla[edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 1]
peakzilla[!edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 2]

peakzilla[edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "0min"]
peakzilla[!edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "20min"]

peakzilla[, Name := NULL]
peakzilla[, edsn := NULL]

rm(dr,files)

######################################################################################

dr <- "inst/mace_analysis"
files <- list.files(dr)
files <- files[grep("bed",files)]

files <- files[grep("border_pair.bed",files,fixed = TRUE)]

mace <- lapply(file.path(dr,files),read.table)
mace <- lapply(mace, data.table)

mace <- lapply(mace,function(x){
  setnames(x,names(x),c("chrom","st","en","border","pval"))
  x[,Name := paste0(chrom,":",st,"-",en)]
  x[,site := mid(IRanges(start = st , end = en))]
  return(x)})

mace <- mapply(function(mosaics,mace){
  idx <- countOverlaps(IRanges(start = mace[,(st)],end = mace[,(en)]),mosaics) > 0
  out <- mace[idx]
  return(out)},mosaics_peaks,mace,SIMPLIFY = FALSE)

mace <- mapply(function(mace,sites){
  mranges <- IRanges(start = mace[,(st)],end = mace[,(en)])
#  ov <- findOverlaps(mranges,sites)  
  reso <- mace[,.(reso = min(abs(site - start(sites)))),by = Name]
  setkey(reso,Name)
#  reso <- reso[as.character(mace[queryHits(ov),(Name)])]
  return(reso)
},mace,sites_mosaics,SIMPLIFY = FALSE)

mace <- mapply(function(x,y)copy(x)[,edsn := y],mace,files,SIMPLIFY = FALSE)
mace <- do.call(rbind,mace)

mace[,edsn := gsub(".border_pair.bed","",edsn)]
mace[,edsn := gsub("Sig70_","",edsn)]

mace[, repl := 0]
mace[, rif := ""]

mace[grep("rep1",edsn),repl := 1]
mace[grep("rep2",edsn),repl := 2]
mace[grep("rif0",edsn),rif := "0min"]
mace[grep("rif20",edsn),rif := "20min"]

mace[,Name := NULL]
mace[,edsn := NULL]

rm(dr, files)

######################################################################################

## gem analysis

dr <- "inst/gem_analysis"
files <- list.files(dr,recursive = TRUE)
files <- files[grep("peaks",files)]
files <- files[grep("bed",files)]
files <- files[grep("insig",files,invert = TRUE)]
files <- files[grep("event",files)]

gem <- lapply(file.path(dr,files),read.table,skip = 1)
gem <- lapply(gem,data.table)

gem <- lapply(gem,function(x){
  setnames(x,names(x),c("chr","st","en","name","signal"))
  x[,site := mid(IRanges(start = st , end = en))]
  return(x)})

gem <- mapply(function(mosaics,gem){
  idx <- countOverlaps(IRanges(start = gem[,(st)],end = gem[,(en)]),mosaics) > 0
  out <- gem[idx]
  return(out)},mosaics_peaks,gem,SIMPLIFY = FALSE)

gem <- mapply(function(gem,sites){
  gemranges <- IRanges(start = gem[,(st)],end = gem[,(en)])
#  ov <- findOverlaps(gemranges,sites)
  reso <- gem[,.(reso = min(abs(site - start(sites)))),by = name]
  setkey(reso,name)
#  reso <- reso[as.character(gem[queryHits(ov),(name)])]
  return(reso)},gem,sites_mosaics,SIMPLIFY = FALSE)
  
gem <- mapply(function(x,y)copy(x)[,edsn := y],gem,basename(files),SIMPLIFY = FALSE)
gem <- do.call(rbind,gem)

gem[,edsn := gsub("_Sig70_with_peaksFDR5_1_GEM_events.bed","",edsn)]
gem[, repl := 0]
gem[, rif := ""]

gem[edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 1]
gem[!edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 2]

gem[edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "0min"]
gem[!edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "20min"]

gem[,name := NULL]
gem[,edsn := NULL]

rm(dr, files)

######################################################################################

dr <- file.path(base_dir,"ChIPexo")
files <- list.files(dr,recursive = TRUE)
files <- files[grep("binding",files)]
files <- files[grep("FDR5",files)]

dpeak <- lapply(file.path(dr,files),read.table,skip = 1)
dpeak <- lapply(dpeak,data.table)

dpeak <- lapply(dpeak,function(x){
  setnames(x,names(x),c("chrID","st","en","peak","strength"))  
  x[,site := mid(IRanges(start = st , end = en))]
  x[,name := paste0(st,"-",en)]
  return(x)})

dpeak <- mapply(function(mosaics,dpeak){
  idx <- countOverlaps(IRanges(start = dpeak[,(st)],end = dpeak[,(en)]),mosaics) > 0
  out <- dpeak[idx]
  return(out)},mosaics_peaks,dpeak,SIMPLIFY = FALSE)

dpeak <- mapply(function(dpeak,sites){
  dranges <- IRanges(start = dpeak[,(st)],end = dpeak[,(en)])
#  ov <- findOverlaps(dranges,sites)
  reso <- dpeak[,.(reso = min(abs(site - start(sites)))),by = name]
  setkey(reso,name)
#  reso <- reso[as.character(dpeak[queryHits(ov),(name)])]
  return(reso)},dpeak,sites_mosaics,SIMPLIFY = FALSE)

    
dpeak <- mapply(function(x,y)copy(x)[,edsn := y],dpeak,basename(files),SIMPLIFY = FALSE)
dpeak <- do.call(rbind,dpeak)

dpeak[,edsn := gsub("_Sig70_sites_G5.txt","",edsn)]
dpeak[, repl := 0]
dpeak[, rif := ""]

dpeak[edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 1]
dpeak[!edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 2]

dpeak[edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "0min"]
dpeak[!edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "20min"]

dpeak[,name := NULL]
dpeak[,edsn := NULL]

rm(dr, files)

resol <- list("peakzilla" = peakzilla,"mace" = mace, "gem" = gem, "dpeak" = dpeak)
resol <- mapply(function(x,y)x[,method := y],resol,names(resol),SIMPLIFY = FALSE)
resol <- do.call(rbind,resol)

pdf(file = file.path(figs_dir,"methods_comparison_resolution_noOv.pdf"))
ggplot(resol, aes(method , reso,fill = method))+geom_boxplot()+facet_grid( repl ~ rif)+
  scale_fill_brewer(palette = "Pastel1")+theme_bw()+theme(legend.position = "none")+
  xlab("Method")+ylab("Resolution")+ylim(0,75)
dev.off()
