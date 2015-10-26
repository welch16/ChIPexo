
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)

######################################################################################

## Condition table

edsn_tab <- function(what){
  stopifnot(what %in% c("exo","pet","set"))
  if(what == "exo"){
    edsn <- as.character(931:938)
    ip <- rep(c("Sig70","SigmaS"),4)
    condition <- rep(c("exp","stat"),each = 4)
    repl <- rep( rep(1:2,each =2),2)
    dt1 <- data.table(edsn,ip,condition,repl)
    edsn <- as.character(1310:1321)
    ip <- rep(c("Beta","Sig70","BetaPF"),4)
    condition <- rep( rep(c("rif0min","rif20min"),each = 3),2)
    repl <- rep(rep(1:2,each = 6))
    dt2 <- data.table(edsn,ip,condition,repl)
    dt <- rbind(dt1,dt2)    
  }else{
    edsn <- as.character(1396:1403)
    ip <- rep(c("Sig70","BetaPF"),4)
    condition <- rep(c("rif0min","rif20min"),each = 2, 2)
    repl <- rep(1:2,each = 4)
    dt <- data.table(edsn,ip,condition,repl)      
  }
  return(dt)
}

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab)
names(char) <- what

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif20min"
bs <- 150
fl <- 150
mc <- detectCores()
figs_dir <- "figs/condition"
fdr <- .01
G <- 1

char <- lapply(char,function(x)x[ip ==tf & condition == rif])

base_dir <- "/p/keles/ChIPexo/volume6/condition"
folder <- paste(tf,rif,sep = "_")

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

figs_dir <- file.path(figs_dir,folder,paste0("FDR",fdr*100))
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

######################################################################################

## load peaks

load_peaks <- function(dir,what,fdr,char)
{
  dir <- file.path(dir,"peaks",what,paste0("FDR",fdr*100))
  files <- list.files(dir)
  edsn <- char[[what]][,(edsn)]
  peaks <- sapply(edsn,function(x)files[grep(x,files)])
  peaks <- file.path(dir,peaks)
  peaks <- lapply(peaks,read.table)
  peaks <- lapply(peaks,data.table)
  peaks <- lapply(peaks,function(x){
    x[,peakId := paste0(V1,":",V2,"-",V3)]
    return(x)})
  names(peaks) <- edsn
  return(peaks)
}

peaks <- lapply(what,function(x)load_peaks(base_dir,x,fdr,char))
names(peaks) <- what

######################################################################################

## load binding sites

load_binding <- function(dir,what,fdr,G,char)
{
  browser()
  dir <- file.path(dir,"binding",what,paste0("FDR",fdr*100),paste0("G_",G))
  files <- list.files(dir)
  edsn <- char[[what]][,(edsn)]
  bs <- sapply(edsn,function(x)files[grep(x,files)])
  bs <- file.path(dir,bs)  
  bs <- lapply(bs,read.table,skip = 1)
  bs <- lapply(bs,data.table)
  bs <- lapply(bs,function(x){
    setnames(x,names(x),c("seqnames","start","end","peakId","strength"))
    return(x)})
  names(bs) <- edsn
  return(bs)
}

bs <- lapply(what,function(x)load_binding(base_dir,x,fdr,G,char))
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
  out <- mapply(merge_peak_sites,peaks[[x]],bs[[x]],MoreArgs = list(sites,x),SIMPLIFY = FALSE)
  return(out)})
names(results) <- what

######################################################################################

res1 <- function(samp,what)
{
  res <- samp[,min(dist),by = peakId]
  setnames(res,names(res),c("peakId","reso"))
  res[,seq := what]
  return(res)
}


resolution <- function(exo,pet,set,repl)
{
  res <- rbind(res1(exo,"exo"),res1(pet,"pet"),res1(set,"set"))
  res[,rep := repl]
  return(res)
}

reso <- mapply(resolution,
              results[["exo"]],results[["pet"]],results[["set"]],1:2,SIMPLIFY = FALSE)
reso <- do.call(rbind,reso)


## resolution plots

pdf(file = file.path(figs_dir,"resolution.pdf"))
ggplot(reso[ rep == 1],aes(seq, reso,colour = seq))+geom_boxplot()+ylim(0,100)+ggtitle("rep1")+
  scale_color_brewer(palette = "Set1")+theme(legend.position = "none")+
  ylab("Resolution")+xlab("")
ggplot(reso[ rep == 2],aes(seq, reso,color = seq))+geom_boxplot()+ylim(0,100)+ggtitle("rep2")+
  scale_color_brewer(palette = "Set1")+theme(legend.position = "none")+
  ylab("Resolution")+xlab("")  
ggplot(reso,aes(seq, reso,colour = seq))+geom_boxplot()+ylim(0,100)+ggtitle("pooled")+
  scale_color_brewer(palette = "Set1")+theme(legend.position = "none")+
  ylab("Resolution")+xlab("")
dev.off()

######################################################################################

sens1 <- function(samp,what,ext)
{
  ## mean different between detected sites
  nPredictions <- samp[,length(unique(site)),by  = peakId]
  nAnnot <- samp[,length(unique(site_start)), by = peakId]
  aveDist <- samp[,mean(diff(unique(site_start))),by = peakId]
  dist2Site <- samp[ ,min(abs(site - site_start )), by  = .(peakId,site_start)]
  idenBS <- dist2Site[, sum(V1 <= ext) , by = peakId]

  setnames(nPredictions,names(nPredictions),c("peakId","nPred"))
  setnames(nAnnot,names(nAnnot),c("peakId","nAnnot"))
  setnames(aveDist,names(aveDist),c("peakId","aveDist"))
  setnames(idenBS,names(idenBS),c("peakId","nIden"))

  out <- merge(nPredictions,nAnnot,by = "peakId")
  out <- merge(out,aveDist,by = "peakId")
  out <- merge(out,idenBS,by = "peakId")

  out[,frac := nIden / nAnnot]

  out[,seq := what]
  return(out)
}

sensitivity <- function(exo,pet,set,repl,ext)
{
  sens <- rbind(sens1(exo,"exo",ext),sens1(pet,"pet",ext),sens1(set,"set",ext))
  sens[,rep := repl]
  return(sens)
}

sens <- mapply(sensitivity ,
           results[["exo"]],
           results[["pet"]],
           results[["set"]],1:2,MoreArgs = list(ext = 20),SIMPLIFY = FALSE)
sens <- do.call(rbind,sens)

xl <- 500
pdf(file = file.path(figs_dir,"sensitivity.pdf"))
ggplot(sens[nAnnot > 1 & between(aveDist,0, xl,incbounds = FALSE) ],
       aes(aveDist,frac,colour = seq,shape = seq))+geom_point()+
  geom_smooth(method = "rlm",se =FALSE,size = 1)+ scale_color_brewer(palette = "Set1")+
  facet_grid(rep ~.) + theme(legend.position = "top")+ggtitle("by replicate")+
  xlim(0,xl)+ylim(0,1)+
  xlab("Average distance between annotated binding sites")+
  ylab("Sensitivity")
dev.off()

