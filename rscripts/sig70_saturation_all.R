
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)

######################################################################################

## Condition table

source("R/base_edsn.R")

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab)
names(char) <- what

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "aerobic"
mc <- 8
figs_dir <- "figs/saturation"
ext <- 20

char <- lapply(char,function(x)x[ip ==tf & condition == rif])

base_dir <- "/p/keles/ChIPexo/volume6/saturation"
folder <- paste(tf,rif,sep = "_")

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

base_dir <- file.path(base_dir,folder)

seeds <- list.files(base_dir)

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

sites[,seqnames := "U00096"]
sites[,strand := ifelse(strand == "F","+","-")]

######################################################################################

load_all <- function(dir,sites)
{
  all_files <- list.files(dir,recursive = TRUE)

  peak_files <- all_files[grep("peaks",all_files)]
  peaks <- lapply(file.path(dir,peak_files),read.table)
  peaks <- lapply(peaks,data.table)
  peaks <- lapply(peaks,function(x){
    x[,peakId := paste0(V1,":",V2,"-",V3)]
    return(x)})

  bs_files <- all_files[grep("binding",all_files)]
  binding <- lapply(file.path(dir,bs_files),read.table,skip = 1)
  binding <- lapply(binding,data.table)
  nms <- c("seqnames","start","end","peakId","strength")
  binding <- lapply(binding,function(x){
    setnames(x,names(x),nms)
    return(x)})

  matchs <- lapply(peak_files,match_peak_sites,bs_files)

  what <- sapply(strsplit(peak_files,"/",fixed = TRUE),function(x)x[2])

  edsn <- sapply(strsplit(peak_files,"/",fixed = TRUE),function(x)x[4])
  sample <- sapply(strsplit(edsn,"_",fixed = TRUE),function(x)x[2])
  sample <- gsub("sample","",sample)
  edsn <- sapply(strsplit(edsn,"_",fixed = TRUE),function(x)x[1])
  
  results <- mcmapply(function(peak,pfile,samp,match,what,bind,sites){
    summaries <- lapply(bind[match],function(x){
      merge_peak_sites(peak,x,what,sites)})
    out <- do.call(rbind,summaries)
    out[, edsn := pfile]
    out[, sample := samp]
    return(out)},peaks,edsn,sample,matchs,what,MoreArgs = list(binding,sites),SIMPLIFY = FALSE,mc.cores =12)

  results <- do.call(rbind,results)

  return(results)

}

merge_peak_sites <- function(peaks,bs,what,sites)
{
  bs <- copy(bs)  
  bs <- merge(bs,peaks[,c("V1","V2","V3","peakId"),with = FALSE],by = "peakId")
  setnames(bs,names(bs),c(names(bs)[1:5],"peak_seqnames","peak_start","peak_end"))

  bs[,peak_size := peak_end - peak_start + 1]

  gr <- bs[,6:8,with = FALSE]
  setnames(gr,names(gr),gsub("peak_","",names(gr)))
  gr <- ChIPUtils::dt2gr(gr)

  ov <- findOverlaps(gr,ChIPUtils::dt2gr(sites))

  sites <- copy(sites[subjectHits(ov)])
  bs <- copy(bs[queryHits(ov)])

  setnames(sites,names(sites),paste0("site_",names(sites)))

  bs <- cbind(bs,sites)
  bs <- bs[order(start)]

  bs[ , site := .5 * (start  + end - 1)]
  bs[, dist := abs(site - site_start )]

  bs[,what := what]

  return(bs)
}

match_peak_sites <- function(peak_file,bs_files)
{
  char <- strsplit(peak_file,"/",fixed = TRUE)[[1]]
  seq <- char[2]
  fdr <- char[3]
  file <- strsplit(char[4],"_")[[1]]
  edsn <- file[1]
  sample <- file[2]
  

  id1 <- grepl(seq,bs_files)
  id2 <- grepl(fdr,bs_files,fixed = TRUE,useBytes = TRUE)  
  ## if(fdr == "FDR1"){
  ##   id2 <- id2 & !grepl("FDR10",bs_files)
  ## }
  id3 <- grepl(edsn,bs_files)
  id4 <- grepl(sample,bs_files)
  return(which( id1 & id4 & id3 ))
}

dirs <- file.path(base_dir,seeds)
summaries <- mclapply(dirs,load_all,sites,mc.cores = 10)

summaries <- mapply(function(x,y)x[,seed := y],summaries,seeds,SIMPLIFY = FALSE)
summaries <- do.call(rbind,summaries)

######################################################################################

npeaks <- summaries[,length(unique(peakId)), by = .(seed,what,sample)]
npeaks[ ,sample := as.numeric(gsub("K","",sample)) ]
p1 <- ggplot(npeaks,aes(sample,V1,colour = what))+geom_point(size = 1)+
  geom_smooth(method = "loess",se = FALSE,size = 1.5)+
  ylim(150,375)+scale_color_brewer(palette = "Set1",name = "Protocol")+
  theme_bw()+theme(legend.position = "top")+
  xlab("Nr. reads (in thousands)")+ylab("Nr. of candidate peaks")

npred <- summaries[,length(peakId),by = .(seed,what,sample)]
npred[ ,sample := as.numeric(gsub("K","",sample)) ]
p2 <- ggplot(npred,aes(sample,V1,colour = what))+geom_point(size = 1)+
  geom_smooth(method = "loess",se = FALSE,size = 1.5)+
  scale_color_brewer(palette = "Set1",name = "Protocol")+
  theme_bw()+theme(legend.position = "top")+ylim(300,1100)+
  xlab("Nr. reads (in thousands)")+ylab("Nr. of predicted events")

ext <- 15
niden <- summaries[ ,sum(dist <= ext),by = .(seed,what,sample)]
niden[ ,sample := as.numeric(gsub("K","",sample)) ]
p3 <- ggplot(niden,aes(sample,V1,colour = what))+geom_point(size = 1)+
  geom_smooth(method = "loess",se = FALSE,size = 1.5)+
  scale_color_brewer(palette = "Set1",name = "Protocol")+
  theme_bw()+theme(legend.position = "top")+
  xlab("Nr. reads (in thousands)")+ylab("Nr. of identificed targets")

reso <- summaries[, min(dist), by = .(seed,what,sample,peakId)]
reso <- reso[,median(V1), by = .(seed,what,sample)]
reso[ ,sample := as.numeric(gsub("K","",sample)) ]

p4 <- ggplot(reso,aes(sample,V1,colour = what))+geom_point(size = 1)+
  geom_smooth(method = "loess",se = FALSE,size = 1.5)+
  scale_color_brewer(palette = "Set1",name = "Protocol")+
  theme_bw()+theme(legend.position = "top")+
  xlab("Nr. reads (in thousands)")+ylab("Resolution")


pdf(file = "figs/for_paper/Sig70_aerobic_npeaks.pdf")
print(p1)
dev.off()

pdf(file = "figs/for_paper/Sig70_aerobic_npredicted_even.pdf")
print(p2)
dev.off()

pdf(file = "figs/for_paper/Sig70_aerobic_nIden_targets.pdf")
print(p3)
dev.off()

pdf(file = "figs/for_paper/Sig70_aerobic_resolution.pdf")
print(p4)
dev.off()

pdf(file = "figs/for_paper/Sig70_aerobic_saturation.pdf")
grid.arrange(p1,p2,p3,p4,nrow = 2)
dev.off()
