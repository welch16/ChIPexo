
rm(list = ls())

library(parallel)
library(data.table)
library(dpeak)
library(GenomicAlignments)

mc <- 24

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
sites[,site := ifelse(str == "F",V2,V3)]

## get peaks
base_dir <- "/p/keles/ChIPexo/volume6/results"
peaks_dir <- file.path(base_dir,"mosaics_peaks","Landick")

get_peaks <- function(peaks_dir,what)
{
  if(what == "exo"){
    ext <- "ChIPexo"
  }else if(what == "pet"){
    ext <- "ChIPseq_PET"
  }else{
    ext <- "ChIPseq_SET"
  }
  dr <- file.path(peaks_dir,ext,"peaks")
  files <- list.files(dr)
  out <- lapply(file.path(dr,files),read.table)
  out <- lapply(out,data.table)
  names(out) <- files
  return(out)
}

exo_peaks <- get_peaks(peaks_dir,"exo")
pet_peaks <- get_peaks(peaks_dir,"pet")
set_peaks <- get_peaks(peaks_dir,"set")


## get binding sites
bs_dir <- file.path(base_dir,"dpeak","Landick")

get_binding_sites <- function(bs_dir,what)
{
  if(what == "exo"){
    ext <- "ChIPexo"
  }else if(what == "pet"){
    ext <- "ChIPseq_PET"
  }else{
    ext <- "ChIPseq_SET"
  }
  dr <- file.path(bs_dir,ext)
  files <- list.files(dr)
  sites <- mclapply(file.path(dr,files),read.table,skip = 1,stringsAsFactors = FALSE,mc.cores = mc)
  sites <- lapply(sites,data.table)
  sites <- lapply(sites,function(x){
    x[,site := floor(.5 * (V2 + V3 - 1))]
    return(x)})
  names(sites) <- files
  return(sites)
}

exo_sites <- get_binding_sites(bs_dir,"exo")
pet_sites <- get_binding_sites(bs_dir,"pet")
set_sites <- get_binding_sites(bs_dir,"set")

## resolution is defined as the distance between regulonDB annotation and it's closest prediction
dt2ir <- function(dt)IRanges(start = dt[[2]],end = dt[[3]])

resolution <- function(dpname,peaks,dpeak_sites,sites,what)
{
  dp <- strsplit(dpname,"_")[[1]][1]

  ## which binding_sites and peaks are we gonna use
  i <- grep(dp,names(peaks))
  j <- grep(dp,names(dpeak_sites))

  ## find peaks that overlap sig70 sites
  site_ov <- findOverlaps(dt2ir(sites),dt2ir(peaks[[i]]))
  sites <- sites[queryHits(site_ov)]
  peak <- peaks[[i]][subjectHits(site_ov)]

  dpeak <- lapply( j ,function(x){
    ov <- findOverlaps(dt2ir(dpeak_sites[[x]]),dt2ir(peak))
    return(dpeak_sites[[x]][queryHits(ov)])})
  names(dpeak) <- names(dpeak_sites[j])

  res <- lapply(dpeak,function(x,sites){
      pos <- x[,(site)]
      res <- sites[,min(abs(site - pos)),by = id]
      setnames(res,names(res),c("id","res"))
      return(res)},sites)

  nms <- names(dpeak)
  comp <- sapply(strsplit(nms,"maxComp"),function(x)x[2])
  comp <- gsub(".bed","",comp)
  comp <- as.integer(comp)

  res <- mapply(function(x,y){
    x[,g := y]
    return(x)},res,comp,SIMPLIFY = FALSE)
  res <- do.call(rbind,res)
  res[,seq := what]

  res[,edsn := gsub("edsn","",dp)]
  
  return(res)
}

exo_res <- mclapply(names(exo_peaks),resolution,exo_peaks,exo_sites,sites,"exo",mc.cores = mc)
pet_res <- mclapply(names(pet_peaks),resolution,pet_peaks,pet_sites,sites,"pet",mc.cores = mc)
set_res <- mclapply(names(set_peaks),resolution,set_peaks,set_sites,sites,"set",mc.cores = mc)

exo_res <- do.call(rbind,exo_res)
pet_res <- do.call(rbind,pet_res)
set_res <- do.call(rbind,set_res)


library(ggplot2)
library(RColorBrewer)
library(grDevices)

figs_dir <- "figs/resolution"

pdf(file = file.path(figs_dir,"resolution_chip_exo_landick.pdf"),width = 12)
ggplot(exo_res,aes(as.factor(g) , res))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+facet_grid(.~edsn )+
  ylab("resolution")+xlab("g")+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_chip_seq_set_landick.pdf"),width = 12)
ggplot(set_res,aes(as.factor(g) , res))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+facet_grid(.~edsn )+
  ylab("resolution")+xlab("g")+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_chip_seq_pet_landick.pdf"),width = 12)
ggplot(pet_res,aes(as.factor(g) , res))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))+facet_grid(.~edsn )+
  ylab("resolution")+xlab("g")+ylim(0,225)
dev.off()

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
    repl <- rep(rep(1:2,each = 6),2)
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

exo_res <- merge(exo_res,exo_char,by = "edsn",allow.cartesian = TRUE)
pet_res <- merge(pet_res,pet_char,by = "edsn",allow.cartesian = TRUE)
set_res <- merge(set_res,set_char,by = "edsn",allow.cartesian = TRUE)

res <- do.call(rbind,list(exo_res,pet_res,set_res))

### new data analysis

pdf(file = file.path(figs_dir,"resolution_Sig70_0min.pdf"))
ggplot(res[res <= 1000 & condition == "rif0min" & ip == "Sig70"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( repl ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_Sig70_20min.pdf"))
ggplot(res[res <= 1000 & condition == "rif20min" & ip == "Sig70"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( repl ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_BetaPF_0min.pdf"))
ggplot(res[res <= 1000 & condition == "rif0min" & ip == "BetaPF"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( repl ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_BetaPF_20min.pdf"))
ggplot(res[res <= 1000 & condition == "rif20min" & ip == "BetaPF"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( repl ~ g)+ylim(0,225)
dev.off()


pdf(file = file.path(figs_dir,"resolution_Sig70_0min_pool.pdf"))
ggplot(res[res <= 1000 & condition == "rif0min" & ip == "Sig70"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( . ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_Sig70_20min_pool.pdf"))
ggplot(res[res <= 1000 & condition == "rif20min" & ip == "Sig70"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( . ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_BetaPF_0min_pool.pdf"))
ggplot(res[res <= 1000 & condition == "rif0min" & ip == "BetaPF"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( . ~ g)+ylim(0,225)
dev.off()

pdf(file = file.path(figs_dir,"resolution_BetaPF_20min_pool.pdf"))
ggplot(res[res <= 1000 & condition == "rif20min" & ip == "BetaPF"],
  aes(seq , res))+geom_boxplot()+
  facet_grid( . ~ g)+ylim(0,225)
dev.off()
