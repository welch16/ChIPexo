
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)

######################################################################################

## Condition table

source("R/base_edsn.R")

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab)
names(char) <- what

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif20min"
bs <- 150
fl <- 150
mc <- 8
figs_dir <- "figs/saturation"
fdr <- .001
G <- 5
k <- 1
ext <- 20

char <- lapply(char,function(x)x[ip ==tf & condition == rif])

base_dir <- "/p/keles/ChIPexo/volume6/saturation"
folder <- paste(tf,rif,sep = "_")

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

figs_dir <- file.path(figs_dir,folder)
check_create(figs_dir)

figs_dir <- file.path(figs_dir,paste0("seed",k))
check_create(figs_dir)

figs_dir <- file.path(figs_dir,paste0("FDR",fdr*100))
check_create(figs_dir)

figs_dir <- file.path(figs_dir,paste0("G_",G))
check_create(figs_dir)

base_dir <- file.path(base_dir,folder)

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

## sites[,start := start - 1]

######################################################################################

## load peaks

extract_nreads <- function(file)
{
  spl <- strsplit(file,"_")[[1]][2]
  spl <- gsub("sample","",spl)
  spl <- gsub(".bed","",spl)
  return(spl)
}


load_peaks <- function(dir,what,seed,fdr,char)
{
  dir <- file.path(dir,paste0("seed",seed),"peaks",what,paste0("FDR",fdr*100))
  files <- list.files(dir)
  edsn <- char[[what]][,(edsn)]
  peaks <- sapply(edsn,function(x)files[grep(x,files)])
  nreads <- sapply(peaks,extract_nreads)
  peaks <- file.path(dir,peaks)
  peaks <- lapply(peaks,read.table)
  peaks <- lapply(peaks,data.table)
  peaks <- lapply(peaks,function(x){
    x[,peakId := paste0(V1,":",V2,"-",V3)]
    return(x)})
  nms <- paste(rep(edsn,each = length(files) / 2), nreads, sep = "_")
  names(peaks) <- nms
  return(peaks)
}

peaks <- lapply(what,function(x)load_peaks(base_dir,x,k,fdr,char))
names(peaks) <- what

######################################################################################

## load binding sites

load_binding <- function(dir,what,seed,fdr,G,char)
{
  # paste0("FDR",fdr*100)
  dir <- file.path(dir,paste0("seed",seed),"binding",what,paste0("G_",G))
  files <- list.files(dir)
  edsn <- char[[what]][,(edsn)]
  bs <- sapply(edsn,function(x)files[grep(x,files)])
  nreads <- sapply(bs,extract_nreads)
  bs <- file.path(dir,bs)  
  bs <- lapply(bs,read.table,skip = 1)
  bs <- lapply(bs,data.table)
  bs <- lapply(bs,function(x){
    setnames(x,names(x),c("seqnames","start","end","peakId","strength"))
    return(x)})
  names(bs) <- paste(rep(edsn,each = length(files) / 2),nreads,sep = "_")
  return(bs)
}

bs <- lapply(what,function(x)load_binding(base_dir,x,k,fdr,G,char))
names(bs) <- what

######################################################################################

## merge peaks, binding events and gold standard sites

merge_peak_sites <- function(peaks,bs,sites,what)
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

  return(bs)
}


results <- lapply(what,function(x){
  out <- mcmapply(merge_peak_sites,peaks[[x]],bs[[x]],MoreArgs = list(sites,x),SIMPLIFY = FALSE,mc.cores = mc)
  return(out)})
names(results) <- what

######################################################################################

## Number of peaks

nr_peaks_by_sample <- function(results)
{
  what <- names(results)
  nr_peaks <- lapply(results,function(x){
    n_peaks <- sapply(x , function(y)y[,length(unique(peakId))])
    return(n_peaks)})
  nr_peaks <- lapply(nr_peaks,function(x){
    nms <- names(x)
    edsn <- sapply(strsplit(nms,"_"),function(z)z[1])
    size <- sapply(strsplit(nms,"_"),function(z)as.numeric(gsub("K","",z[2]))*1e3)
    out <- data.table(edsn = factor(edsn),size,nr_peaks = x,repl = 0)
    return(out)})
  nr_peaks <- mapply(function(x,y)x[,seq := y],nr_peaks,what,SIMPLIFY = FALSE)
  nr_peaks <- lapply(nr_peaks,function(x)x[,repl := as.numeric(edsn)])
  nr_peaks <- do.call(rbind,nr_peaks)  
  return(nr_peaks)
}

nr_cand_region <- nr_peaks_by_sample(results)

pdf(file = file.path(figs_dir,"Nr_candidate_regions.pdf"),width = 6,height = 5)
ggplot(nr_cand_region[,mean(nr_peaks),by = .(size,seq)],aes(size , V1,colour = seq))+
  geom_line(size = 1)+theme(legend.position = "top")+ylim(140,400)+
  scale_color_brewer(name = "",palette = "Set1")+ylab("Number of candidate regions")+xlab("Number of reads")
ggplot(nr_cand_region,aes(size , nr_peaks,colour = seq,linetype = as.factor(repl)))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1",guide = "none")+scale_linetype_manual(values = c(1,2),name = "rep")+
  ylab("Number of candidate regions")+xlab("Number of reads")+theme(legend.position = "top")+ylim(140,400)
dev.off()

######################################################################################

## Number of predicted events

nr_predicted_events_by_sample <- function(results)
{
  what <- names(results)
  nr_events <- lapply(results,function(x){
    n_ev <- sapply(x,nrow)
    return(n_ev)})
  
  nr_events <- lapply(nr_events,function(x){
    nms <- names(x)
    edsn <- sapply(strsplit(nms,"_"),function(z)z[1])
    size <- sapply(strsplit(nms,"_"),function(z)as.numeric(gsub("K","",z[2]))*1e3)
    out <- data.table(edsn = factor(edsn),size,nr_events = x,repl = 0)
    return(out)})
  nr_events <- mapply(function(x,y)x[,seq := y],nr_events,what,SIMPLIFY = FALSE)
  nr_events <- lapply(nr_events,function(x)x[,repl := as.numeric(edsn)])
  nr_events <- do.call(rbind,nr_events)  
  return(nr_events)
}

nr_pred_events <- nr_predicted_events_by_sample(results)

pdf(file = file.path(figs_dir,"Nr_predicted_events.pdf"),width = 6,height = 5)
ggplot(nr_pred_events[,mean(nr_events),by = .(size,seq)],aes(size , V1,colour = seq))+
  geom_line(size = 1)+theme(legend.position = "top")+ylim(250,1e3)+
  scale_color_brewer(name = "",palette = "Set1")+ylab("Number of predicted events")+xlab("Number of reads")
ggplot(nr_pred_events,aes(size , nr_events,colour = seq,linetype = as.factor(repl)))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1",guide = "none")+scale_linetype_manual(values = c(1,2),name = "rep")+
  ylab("Number of predicted events")+xlab("Number of reads")+theme(legend.position = "top")+ylim(250,1e3)
dev.off()

######################################################################################

## Number of identified targets

nr_iden_target_by_sample <- function(results,ext)
{
  what <- names(results)
  nr_iden <- lapply(results,function(x){
    n_id <- sapply(x,function(y)nrow(y[dist <= ext]))
    return(n_id)})  
  nr_iden <- lapply(nr_iden,function(x){
    nms <- names(x)
    edsn <- sapply(strsplit(nms,"_"),function(z)z[1])
    size <- sapply(strsplit(nms,"_"),function(z)as.numeric(gsub("K","",z[2]))*1e3)
    out <- data.table(edsn = factor(edsn),size,nr_id_target = x,repl = 0)
    return(out)})
  
  nr_iden <- mapply(function(x,y)x[,seq := y],nr_iden,what,SIMPLIFY = FALSE)
  nr_iden <- lapply(nr_iden,function(x)x[,repl := as.numeric(edsn)])
  nr_iden <- do.call(rbind,nr_iden)

  return(nr_iden)

}

nr_iden_targets <- nr_iden_target_by_sample(results,ext)

pdf(file = file.path(figs_dir,"Nr_iden_targets.pdf"),width = 6,height = 5)
ggplot(nr_iden_targets[,mean(nr_id_target),by = .(size,seq)],aes(size , V1,colour = seq))+
  geom_line(size = 1)+theme(legend.position = "top")+ylim(50,350)+
  scale_color_brewer(name = "",palette = "Set1")+ylab("Number of identified targets")+xlab("Number of reads")
ggplot(nr_iden_targets,aes(size , nr_id_target,colour = seq,linetype = as.factor(repl)))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1",guide = "none")+scale_linetype_manual(values = c(1,2),name = "rep")+
  ylab("Number of identified targets")+xlab("Number of reads")+theme(legend.position = "top")+ylim(50,350)
dev.off()


######################################################################################

## Resolution

resolution_by_sample <- function(results)
{
  what <- names(results)
  reso <- lapply(results,function(x){
    out <- sapply(x,function(y){
      dist <- y[, min(dist), by = peakId][,(V1)]
      return(median(dist))})
    return(out)})
      
  reso <- lapply(reso,function(x){
    nms <- names(x)
    edsn <- sapply(strsplit(nms,"_"),function(z)z[1])
    size <- sapply(strsplit(nms,"_"),function(z)as.numeric(gsub("K","",z[2]))*1e3)
    out <- data.table(edsn = factor(edsn),size,res = x,repl = 0)
    return(out)})
  
  reso <- mapply(function(x,y)x[,seq := y],reso,what,SIMPLIFY = FALSE)
  reso <- lapply(reso,function(x)x[,repl := as.numeric(edsn)])
  reso <- do.call(rbind,reso)

  return(reso)
}

reso <- resolution_by_sample(results)

pdf(file = file.path(figs_dir,"Resolution.pdf"),width = 6,height = 5)
ggplot(reso[,mean(res),by = .(size,seq)],aes(size , V1,colour = seq))+
  geom_line(size = 1)+theme(legend.position = "top")+ylim(0,25)+
  scale_color_brewer(name = "",palette = "Set1")+ylab("Resolution")+xlab("Number of reads")
ggplot(reso,aes(size , res,colour = seq,linetype = as.factor(repl)))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1",guide = "none")+scale_linetype_manual(values = c(1,2),name = "rep")+
  ylab("Resolution")+xlab("Number of reads")+theme(legend.position = "top")+ylim(0,25)
dev.off()
