
rm(list = ls())

library(parallel)
library(data.table)
library(dpeak)
library(GenomicAlignments)

mc <- 24
promoters <- "/p/keles/ChIPexo/volume6/gold_standard/Landick/PROMOTER.txt"
dt <- read.delim(promoters,sep= "\t", comment.char = "#",header = FALSE)
dt <- data.table(dt)
dt <- dt[,1:5,with = FALSE]
setnames(dt,names(dt),c("id","name","strand","pos","Sigma"))
sites <- dt[grep("Sigma70",Sigma)]


site_dir_exo <- "/p/keles/ChIPexo/volume6/results/dpeak/Landick/ChIPexo/"
files <- list.files(site_dir_exo)

dpeak_sites <- mclapply(file.path(site_dir_exo,files),read.table ,skip = 1,
  stringsAsFactors = FALSE,mc.cores = mc)                      

dpeak_sites <- lapply(dpeak_sites,data.table)
names(dpeak_sites) <- files

nms <- c("seqnames","siteStart","siteEnd","name","strength")
dpeak_sites <- lapply(dpeak_sites,function(x,nms){
  setnames(x,names(x),nms)
  return(x)},nms)

dpeak_sites <- lapply(dpeak_sites,function(x){
  x[,site := floor(.5 * (siteEnd + siteStart))]
  return(x)})

## resolution is defined as the distance between regulonDB annotation and it's closest prediction

resolution <- function(dpeak,sites)
{
  dpeak_pos <- dpeak[,(site)]
  res <- sites[, min(abs(dpeak_pos - pos)),by = id]
  setnames(res,names(res),c("id","res"))
  return(res)
}


res <- lapply(dpeak_sites,resolution,sites)

library(ggplot2)
library(RColorBrewer)
library(grDevices)

figs_dir <- "figs/resolution"

res <- mapply(function(x,y){
  x[,sample := y]
  return(x)},res,names(res),SIMPLIFY = FALSE)

res <- do.call(rbind,res)

pdf(file = file.path(figs_dir,"resolution_chip_exo_landick.pdf"))
ggplot(res,aes(sample , res))+geom_boxplot()+ylim(0,150)+
  theme(axis.text.x = element_text(angle = 90))
dev.off()





