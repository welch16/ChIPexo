
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)

######################################################################################

## Initial parameters

mc <- detectCores()
figs_dir <- "figs/condition"

base_dir <- "/p/keles/ChIPexo/volume6/condition"

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)


figs_dir <- file.path(figs_dir,"all_together")
check_create(figs_dir)

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

conditions <-list.files(base_dir)
dirs <- file.path(base_dir,conditions)

fdr <- c(.1,1,10,5)
what <- c("exo","pet","set")

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
  edsn <- sapply(strsplit(edsn,"_",fixed = TRUE),function(x)x[1])
  
  results <- mcmapply(function(peak,pfile,match,what,bind,sites){
    summaries <- lapply(bind[match],function(x){
      merge_peak_sites(peak,x,what,sites)})
    summaries <- mapply(function(x,y,z)x[,maxG := y],summaries,c(1,3,5),SIMPLIFY = FALSE)
    out <- do.call(rbind,summaries)
    out[, edsn := pfile]
    return(out)},peaks,edsn,matchs,what,MoreArgs = list(binding,sites),SIMPLIFY = FALSE,mc.cores =12)

  fdr <- sapply(strsplit(peak_files,"/",fixed = TRUE),function(x)x[3])
  results <- mapply(function(x,y)x[,FDR := as.numeric(gsub("FDR","",y))],results,fdr,SIMPLIFY = FALSE)
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
  edsn <- strsplit(char[4],"_")[[1]][1]

  id1 <- grepl(seq,bs_files)
  id2 <- grepl(fdr,bs_files,fixed = TRUE,useBytes = TRUE)
  if(fdr == "FDR1"){
    id2 <- id2 & !grepl("FDR10",bs_files)
  }
  id3 <- grepl(edsn,bs_files)
  return(which( id1 & id2 & id3 ))
}

summaries <- lapply(dirs,load_all,sites)

summaries <- mapply(function(x,y)x[,experiment := y],summaries,conditions,SIMPLIFY = FALSE)
summaries <- do.call(rbind,summaries)

######################################################################################

library(scales)

reso <- summaries[,min(dist),by = .(maxG,FDR,experiment,what,peakId,edsn)]

reso[ , repl := 0]
reso[ experiment == "Sig70_aerobic" & what == "exo", repl := ifelse(edsn == "edsn931",1,2)]
reso[ experiment == "Sig70_aerobic" & what == "pet", repl := ifelse(edsn == "edsn788",1,2)]
reso[ experiment == "Sig70_aerobic" & what == "set", repl := ifelse(edsn == "edsn80",1,2)]
reso[ experiment == "Sig70_rif0min" & what == "exo", repl := ifelse(edsn == "edsn1311",1,2)]
reso[ experiment == "Sig70_rif0min" & what != "exo", repl := ifelse(edsn == "edsn1396",1,2)]
reso[ experiment == "Sig70_rif20min" & what == "exo", repl := ifelse(edsn == "edsn1314",1,2)]
reso[ experiment == "Sig70_rif20min" & what != "exo", repl := ifelse(edsn == "edsn1398",1,2)]





pdf(file = file.path(figs_dir,"Resolution_all_together.pdf"))
ggplot(reso,aes( as.factor(maxG) , V1,fill =as.factor(FDR)))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(experiment ~ what)+scale_fill_brewer(palette = "Pastel1",name = "FDR")+
  theme_bw()+xlab("Max. nr. of binding sites")+ylab("Resolution")+
  theme(legend.position = "top")
dev.off()

pdf(file = file.path(figs_dir,"Resolution_maxG5.pdf"))
ggplot(reso[maxG == 5],aes( what , V1,fill =as.factor(FDR)))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(experiment ~ .)+scale_fill_brewer(palette = "Pastel1",name = "FDR")+
  theme_bw()+xlab("Sequencing protocol")+ylab("Resolution")+
  theme(legend.position = "top")
dev.off()


pdf(file = file.path(figs_dir,"Resolution_FDR5_maxG5.pdf"))
ggplot(reso[maxG == 5 & FDR == 5],aes( what , V1,fill = what))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(experiment ~ .)+scale_fill_brewer(palette = "Pastel2")+
  theme_bw()+xlab("Sequencing protocol")+ylab("Resolution")+
  theme(legend.position = "none")
dev.off()


pdf(file = file.path(figs_dir,"Resolution_FDR5_maxG5_aerobic.pdf"),width =6,height = 4)
ggplot(reso[maxG == 5 & FDR == 5 & experiment == "Sig70_aerobic"],aes( what , V1,fill = what))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(. ~ experiment)+scale_fill_brewer(palette = "Pastel2")+
  theme_bw()+xlab("Sequencing protocol")+ylab("Resolution")+
  theme(legend.position = "none")
dev.off()


pdf(file = file.path(figs_dir,"Resolution_FDR5_maxG5_rif0.pdf"),width =6,height = 4)
ggplot(reso[maxG == 5 & FDR == 5 & experiment == "Sig70_rif0min"],aes( what , V1,fill = what))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(. ~ repl)+scale_fill_brewer(palette = "Pastel2")+
  theme_bw()+xlab("Sequencing protocol")+ylab("Resolution")+
  theme(legend.position = "none")
dev.off()

pdf(file = file.path(figs_dir,"Resolution_FDR5_maxG5_rif20.pdf"),width =6,height = 4)
ggplot(reso[maxG == 5 & FDR == 5 & experiment == "Sig70_rif20min"],aes( what , V1,fill = what))+
  geom_boxplot(outlier.size = .8,lwd = .2)+ylim(0,35)+
  facet_grid(. ~ repl)+scale_fill_brewer(palette = "Pastel2")+
  theme_bw()+xlab("Sequencing protocol")+ylab("Resolution")+
  theme(legend.position = "none")
dev.off()


######################################################################################

sensitivity <- function(summaries,condition,G,fdr,minStrength,ext)
{
  data <- copy(summaries[experiment == condition & maxG == G])
  data <- split(data,data$what)

  peaks <- lapply(data,function(x)x[,.(seqnames,peak_start,peak_end)])
  peaks <- lapply(peaks,function(x){
    setnames(x,names(x),gsub("peak_","",names(x)))
    return(x)})
  peaks <- lapply(peaks,dt2gr)
  exo_subset <- which(countOverlaps(peaks[[1]],peaks[[2]]) > 0 & countOverlaps(peaks[[1]],peaks[[3]]) > 0)
  pet_subset <- which(countOverlaps(peaks[[2]],peaks[[1]]) > 0 & countOverlaps(peaks[[2]],peaks[[3]]) > 0)
  set_subset <- which(countOverlaps(peaks[[3]],peaks[[2]]) > 0 & countOverlaps(peaks[[3]],peaks[[1]]) > 0)

  data[["exo"]] <- data[["exo"]][exo_subset][strength > minStrength]
  data[["pet"]] <- data[["pet"]][pet_subset][strength > minStrength]
  data[["set"]] <- data[["set"]][set_subset][strength > minStrength]

  data <- mapply(sens1,data,what,MoreArgs = list(ext = ext),SIMPLIFY =FALSE)
  data <- lapply(data,function(x)x[!is.nan(aveDist)])

  data <- do.call(rbind,data)
  data[,experiment := "condition"]
  data[,maxG := G]
  data[,FDR:=fdr]

  return(data)
}

sens1 <- function(samp,what,ext)
{
  browser()
  ## mean different between detected sites
  nPredictions <- samp[,length(unique(site)),by  = peakId]
  nAnnot <- samp[,length(unique(site_start)), by = peakId]
  aveDist <- samp[,mean(abs(diff(unique(site_start)))),by = peakId]
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

library(ChIPUtils)
library(MASS)

condition <- "Sig70_aerobic"
fdr <- 10
G <- 3
minStrength <- 3e3
ext <- 20
sens <- sensitivity(summaries,condition,G,fdr,minStrength,ext)

ggplot(sens,aes(aveDist,frac,colour = seq))+geom_point()+geom_smooth(method = "rlm",se = FALSE)+
  ylim(0,1)+xlim(0,300)
dev.off()
       
## sensitivity <- function(exo,pet,set,repl,ext)
## {
##   sens <- rbind(sens1(exo,"exo",ext),sens1(pet,"pet",ext),sens1(set,"set",ext))
##   sens[,rep := repl]
##   return(sens)
## }

## sens <- mapply(sensitivity ,
##            results[["exo"]],
##            results[["pet"]],
##            results[["set"]],1:2,MoreArgs = list(ext = 20),SIMPLIFY = FALSE)
## sens <- do.call(rbind,sens)

## xl <- 400
## pdf(file = file.path(figs_dir,"sensitivity.pdf"))
## dev.off()

