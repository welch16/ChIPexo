
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)

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

sites[,seqnames := "U00096"]
sites[,strand := ifelse(strand == "F","+","-")]

######################################################################################

match_peak_sites <- function(peak_file,bs_files)
{
  
  char <- strsplit(peak_file,"/",fixed = TRUE)[[1]]
  fdr <- char[2]
  edsn <- strsplit(char[3],"_")[[1]][1]

  id2 <- grepl(fdr,bs_files,fixed = TRUE,useBytes = TRUE)
  if(fdr == "FDR1"){
    id2 <- id2 & !grepl("FDR10",bs_files)
  }
  id3 <- grepl(edsn,bs_files)
  return(which(  id2 & id3 ))
}

merge_peak_sites <- function(peaks,bs,fdr,sites)
{
  bs <- copy(bs)  
  bs <- merge(bs,peaks[,c("V1","V2","V3","peakId"),with = FALSE],by = "peakId")
  setnames(bs,names(bs),c(names(bs)[1:5],"peak_seqnames","peak_start","peak_end"))

  bs[,peak_size := peak_end - peak_start + 1]

  gr <- bs[,6:8,with = FALSE]
  setnames(gr,names(gr),gsub("peak_","",names(gr)))
  gr <- ChIPUtils::dt2gr(gr)

  ov <- findOverlaps(ranges(gr),ranges(ChIPUtils::dt2gr(sites)))

  sites <- copy(sites[subjectHits(ov)])
  bs <- copy(bs[queryHits(ov)])

  setnames(sites,names(sites),paste0("site_",names(sites)))

  bs <- cbind(bs,sites)
  bs <- bs[order(start)]

  bs[ , site := .5 * (start  + end - 1)]
  bs[, dist := min(abs(site - site_start ),
         abs(site - site_end))]

  bs[,fdr := fdr]

  return(bs)
}

load_all <- function(dir,sites)
{
  
  all_files <- list.files(dir,recursive = TRUE)

  peak_files <- all_files[grep("peaks",all_files)]
  peaks <- lapply(file.path(dir,peak_files),read.table)
  peaks <- lapply(peaks,data.table)
  peaks <- lapply(peaks,function(x){
    x[,peakId := paste0(V1,":",V2,"-",V3)]
    return(x)})

  bs_files <- all_files[grep("sites",all_files)]
  binding <- lapply(file.path(dir,bs_files),read.table,skip = 1)
  binding <- lapply(binding,data.table)
  nms <- c("seqnames","start","end","peakId","strength")
  binding <- lapply(binding,function(x){
    setnames(x,names(x),nms)
    return(x)})

  matchs <- lapply(peak_files,match_peak_sites,bs_files)
  

  fdr <- sapply(strsplit(peak_files,"/",fixed = TRUE),function(x)x[2])

  edsn <- sapply(strsplit(peak_files,"/",fixed = TRUE),function(x)x[3])
  edsn <- sapply(strsplit(edsn,"_",fixed = TRUE),function(x)x[1])
  
  results <- mcmapply(function(peak,pfile,match,what,bind,sites){
    summaries <- lapply(bind[match],function(x){
      merge_peak_sites(peak,x,what,sites)})
##    summaries <- mapply(function(x,y,z)x[,maxG := y],summaries,c(1,3,5),SIMPLIFY = FALSE)   
    out <- do.call(rbind,summaries)
    out[, edsn := pfile]
    return(out)},peaks,edsn,matchs,fdr,MoreArgs = list(binding,sites),SIMPLIFY = FALSE,mc.cores =12)


  results <- do.call(rbind,results)

  return(results)

}


dirs <- file.path(base_dir,list.files(base_dir))
dirs <- dirs[grep("inputs",dirs,invert = TRUE)]

summaries <- lapply(dirs,load_all,sites)

summaries <- mapply(function(x,y)x[,seq := y],summaries,basename(dirs),SIMPLIFY = FALSE)
summaries <- do.call(rbind,summaries)


summaries[, repl := 0]
summaries[, rif := ""]


summaries[edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 1]
summaries[!edsn %in% c("edsn1311","edsn1314","edsn1396","edsn1398"), repl := 2]

summaries[edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "0min"]
summaries[!edsn %in% c("edsn1311","edsn1317","edsn1396","edsn1400"), rif := "20min"]

######################################################################################

### resolution boxplots

library(scales)

reso <- summaries[,.(reso = min(dist)),by = .(fdr,seq,peakId,edsn,rif,edsn,repl)]
reso[,edsn := gsub("edsn","",edsn)]

pdf(file = file.path(figs_dir,"resolution_fdr5_all.pdf"))
p <- ggplot(reso[reso <= 100], aes( edsn , reso,colour = seq))+ geom_boxplot() +
  facet_grid(. ~ seq, scales = "free_x")+
  ylim(0,75)+theme_bw()+theme(axis.text.x = element_text(angle = 90),legend.position = "none")+
  scale_color_brewer(palette = "Set1")+ylab("resolution")
print(p)
dev.off()

pdf(file = file.path(figs_dir,"resolution_fdr5_rif0min.pdf"))
p <- ggplot(reso[rif == "0min" & reso <= 100], aes( as.factor(repl), reso,colour = seq))+
  geom_boxplot() + facet_grid(. ~ seq, scales = "free_x")+
  ylim(0,75)+theme_bw()+theme(legend.position = "none")+
  scale_color_brewer(palette = "Set1")+xlab("replicate")+ylab("resolution")
print(p)
dev.off()

pdf(file = file.path(figs_dir,"resolution_fdr5_rif20min.pdf"))
p <- ggplot(reso[rif == "20min" & reso <= 100], aes( as.factor(repl), reso,colour = seq))+
  geom_boxplot() + facet_grid(. ~ seq, scales = "free_x")+
  ylim(0,75)+theme_bw()+theme(legend.position = "none")+
  scale_color_brewer(palette = "Set1")+xlab("replicate")+ylab("resolution")
print(p)
dev.off()


######################################################################################

### sensitivity analysis
## the idea of this plot is to check that as the distance between binding sites
## increases, then we detect more of the binding events that form a peak.

library(ChIPUtils)
library(MASS)


sensitivity <- function(summaries,FDR,minStrength,ext,rifm = NULL,repli = NULL)
{

  data <- copy(summaries[fdr == FDR])
  if(!is.null(rifm)){
    data <- data[rif == rifm]
  }
  if(!is.null(repli)){
    data <- data[repl == repli]
  }
  
  data <- split(data,data$seq)

  peaks <- lapply(data,function(x)x[,.(seqnames,peak_start,peak_end)])
  peaks <- lapply(peaks,function(x){
    setnames(x,names(x),gsub("peak_","",names(x)))
    return(x)})
  peaks <- lapply(peaks,dt2gr)
  
  exo_subset <- which(countOverlaps(peaks[[1]],peaks[[2]]) > 0 & countOverlaps(peaks[[1]],peaks[[3]]) > 0)
  pet_subset <- which(countOverlaps(peaks[[2]],peaks[[1]]) > 0 & countOverlaps(peaks[[2]],peaks[[3]]) > 0)
  set_subset <- which(countOverlaps(peaks[[3]],peaks[[2]]) > 0 & countOverlaps(peaks[[3]],peaks[[1]]) > 0)

  data[["ChIPexo"]] <- data[["ChIPexo"]][exo_subset][strength > minStrength]
  data[["ChIPseq_PET"]] <- data[["ChIPseq_PET"]][pet_subset][strength > minStrength]
  data[["ChIPseq_SET"]] <- data[["ChIPseq_SET"]][set_subset][strength > minStrength]

  data <- mapply(sens1,data,names(data),MoreArgs = list(ext = ext),SIMPLIFY =FALSE)
  
  data <- lapply(data,function(x)x[!is.nan(aveDist)])
  data <- do.call(rbind,data)  
  data[,fdr := FDR]
  data[,rif := rifm]
  data[,repl := repli]

  return(data)
}


nr_identified <- function(site_start,site_end,pred,ext)
{
  upred <- unique(pred)
  ustart <- unique(site_start)
  uend <- unique(site_end)
  
  M <- matrix(0,  nrow = length(upred),ncol = length(ustart))
  for( i in 1:nrow(M)){
    for(j in 1:ncol(M)){
      M[i,j] <- min(abs(upred[i] - ustart[j]),abs(upred[i] - uend[j]))
    }
  }
  
  a = rep(0,ncol(M))
  for(j in 1:ncol(M))if(any(M[,j] <= ext))a[j]=1

  
  return(sum(a))

}



sens1 <- function(samp,what,ext)
{

  ## mean different between detected sites

  out <- samp[,.(nPred = length(unique(site)),
                 nAnnot = length(unique(site_start)),
                 aveDist = mean(abs(diff(unique(site_start)))),
                 nIden = nr_identified(site_start,site_end,site,ext)
                 ),by = peakId]

  ## setkey(samp,peakId)
  ## setkey(out,peakId)
  ## zz = out[,(peakId)][out[,which(!is.na(aveDist))]]

  ## samp[ zz[3], nr_identified(site_start,site_end,site,ext),by = peakId]
  ## samp[ zz[1], nr_identified(site_start,site_end,site,ext),by = peakId]

  out[,fracAnnot := nIden / nAnnot]
  out[,fracPred := nIden / nPred]
  out[,seq := what]


  ## nPredictions <- samp[,length(unique(site),by  = peakId]
  ## nAnnot <- samp[,length(unique(site_start)), by = peakId]
  ## aveDist <- samp[,mean(abs(diff(unique(site_start)))),by = peakId]
  ## dist2Site <- samp[ ,min(abs(site - site_start )), by  = .(peakId,site_start)]
  ## idenBS <- dist2Site[, sum(V1 <= ext) , by = peakId]

  ## setnames(nPredictions,names(nPredictions),c("peakId","nPred"))
  ## setnames(nAnnot,names(nAnnot),c("peakId","nAnnot"))
  ## setnames(aveDist,names(aveDist),c("peakId","aveDist"))
  ## setnames(idenBS,names(idenBS),c("peakId","nIden"))

  ## out <- merge(nPredictions,nAnnot,by = "peakId")
  ## out <- merge(out,aveDist,by = "peakId")
  ## out <- merge(out,idenBS,by = "peakId")

  return(out)
}

minStrength <- 1e3
ext <- 10
FDR <- "FDR5"
sens <- sensitivity(summaries,FDR,minStrength,ext)
ggplot(sens[!(aveDist > 200 & nIden == 0)],aes(aveDist,fracPred,colour = seq))+geom_point()+
  geom_smooth(method = "rlm",se = FALSE,fulrange = TRUE)+
  ylim(0,1)+xlim(0,400)+scale_color_brewer(palette = "Set1")+theme_bw()+
  theme(legend.position = "top")
dev.off()

## ggplot(summaries,aes(strength,colour = seq))+stat_density(geom = "line")+
##   facet_grid( repl  ~ rif)+
##   scale_x_log10()+scale_color_brewer(palette = "Set1")+
##   theme_bw()+theme(legend.position = "top"); dev.off()



