
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
  
exo_char <- edsn_tab("exo")
pet_char <- edsn_tab("pet")
set_char <- edsn_tab("set")

######################################################################################

## Initial parameters

tf <- "Sig70"
rif <- "rif20min"
bs <- 150
fl <- 150
## fdr <- .25
## thresh <- 10
mc <- detectCores()
## g <- 5
figs_dir <- "figs/condition"

exo <- exo_char[ip == tf & condition == rif ]
pet <- pet_char[ip == tf & condition == rif ]
set <- set_char[ip == tf & condition == rif ]

base_dir <- "/p/keles/ChIPexo/volume6/condition"
folder <- paste(tf,rif,sep = "_")

figs_dir <- file.path(figs_dir,folder)
if(!dir.exists(figs_dir))dir.create(figs_dir)

base_dir <- file.path(base_dir,folder)

if(!dir.exists(base_dir)){
  dir.create(base_dir)
}

what <- c("exo","pet","set")
bases <- lapply(what,function(x)file.path(base_dir,x))

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

load_peaks <- function(dir)
{
  files <- list.files(dir)
  files <- files[grep("peaks",files)]
  peaks <- lapply(file.path(dir,files),read.table)
  peaks <- lapply(peaks,data.table)
  peaks <- lapply(peaks,function(x){
    x[,peakId := paste0(V1,":",V2,"-",V3)]
    return(x)})
  return(peaks)
}

exo_peaks <- load_peaks(bases[[1]])
pet_peaks <- load_peaks(bases[[2]])
set_peaks <- load_peaks(bases[[3]])

######################################################################################

## load binding sites

load_binding <- function(dir)
{
  files <- list.files(dir)
  files <- files[grep("bs",files)]
  bs <- lapply(file.path(dir,files),read.table,skip = 1)
  bs <- lapply(bs,data.table)
  bs <- lapply(bs,function(x){
    setnames(x,names(x),c("seqnames","start","end","peakId","strength"))
    return(x)})
  return(bs)
}

exo_bs <- load_binding(bases[[1]])
pet_bs <- load_binding(bases[[2]])
set_bs <- load_binding(bases[[3]])

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

exo_results <- mapply(merge_peak_sites, exo_peaks,exo_bs,MoreArgs = list(sites , "exo"),SIMPLIFY = FALSE)
pet_results <- mapply(merge_peak_sites, pet_peaks,pet_bs,MoreArgs = list(sites , "pet"),SIMPLIFY = FALSE)
set_results <- mapply(merge_peak_sites, set_peaks,set_bs,MoreArgs = list(sites , "set"),SIMPLIFY = FALSE)

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

reso <- mapply(resolution , exo_results,pet_results,set_results,1:2,SIMPLIFY = FALSE)
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


k <- "U00096:1014450-1015949"
k1 <- "U00096:0-449"
k2 <- "U00096:4637850-4639499"

sensitivity <- function(exo,pet,set,repl,ext)
{
  sens <- rbind(sens1(exo,"exo",ext),sens1(pet,"pet",ext),sens1(set,"set",ext))
  sens[,rep := repl]
  return(sens)
}

sens <- mapply(sensitivity , exo_results,pet_results,set_results,1:2,MoreArgs = list(ext = 20),SIMPLIFY = FALSE)
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

## ggplot(sens[nAnnot > 1 & nIdBS > 1],aes(aveDist,frac,colour = seq))+geom_point()+
##   geom_smooth(method = "rlm",se =FALSE,size = 1)+ scale_color_brewer(palette = "Set1")+
##   theme(legend.position = "top")+ggtitle("pooled")+xlim(0,500)

