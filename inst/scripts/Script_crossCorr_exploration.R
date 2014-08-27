
rm(list = ls())
library(parallel)
library(spp)
library(ChIPQC)
library(reshape2)
library(ggplot2)

# This script only converts the data from bam file format to RData

dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder
ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)
names(ff) = folder

ma <- function(x,n)filter(x,rep(1/n,n),sides = 2)

read_chipqc <- function(files,mc.cores)
{
  out = mclapply(files,function(x)ChIPQCsample(x,annotation = NULL,runCrossCor = TRUE),mc.cores = mc.cores)
  return(out)
}

read_spp <- function(files,mc.cores)
{
  out = mclapply(files,function(x)read.bam.tags(x),mc.cores =8)
  return(out)
}

mc.cores = 8
exo_chipqc = read_chipqc(ff[[1]],mc.cores)
pet_chipqc = read_chipqc(ff[[2]],mc.cores)
set_chipqc = read_chipqc(ff[[3]],mc.cores)

exo_spp = read_spp(ff[[1]],mc.cores)
pet_spp = read_spp(ff[[2]],mc.cores)
set_spp = read_spp(ff[[3]],mc.cores)

cross_corr_spp <- function(seq_spp,range,bin)
{  
  out = lapply(seq_spp,function(x,range,bin)get.binding.characteristics(x,srange= range,bin = bin),range,bin)
  return(out)    
}

range = c(1,300)
bin = 1
exo_crossCorr_spp = cross_corr_spp(exo_spp,range,bin)
pet_crossCorr_spp = cross_corr_spp(pet_spp,range,bin)
set_crossCorr_spp = cross_corr_spp(set_spp,range,bin)

shift = seq(range[1],range[2],by=bin)

cross_corr_df <- function(seq_chipqc,seq_crossCorr_spp,seq_names,shift,mc.cores)
{
  n = length(seq_chipqc) 
  stopifnot( n == length(seq_crossCorr_spp))
  stopifnot( n == length(seq_names))
  df_list = mclapply(1:n,function(i,seq_chipqc,seq_crossCorr_spp,seq_names,shift){
    qc_crossCorr = seq_chipqc[[i]]@CrossCorrelation
    spp_crossCorr = seq_crossCorr_spp[[i]]$cross.correlation$y
    df = melt(data.frame(qc = qc_crossCorr,spp=spp_crossCorr))
    names(df) = c("method","crossCorr")
    df$shift = rep(shift,2)
    df$sample = seq_names[i]
    return(df)
  },seq_chipqc,seq_crossCorr_spp,seq_names,shift,mc.cores =mc.cores)
  return(df_list)
}

exo_df_list = cross_corr_df(exo_chipqc,exo_crossCorr_spp,files[[1]],shift,mc.cores)
pet_df_list = cross_corr_df(pet_chipqc,pet_crossCorr_spp,files[[2]],shift,mc.cores)
set_df_list = cross_corr_df(set_chipqc,set_crossCorr_spp,files[[3]],shift,mc.cores)

save(list = c("exo_df_list","pet_df_list","set_df_list"),file = "data/cross_corr.RData")

## df = do.call(rbind,set_df_list)
## df$sample = factor(df$sample)
## p =ggplot(df,aes(shift,crossCorr,colour = method))+geom_line()+
##   facet_wrap(~sample,scales = "free",nrow = length(levels(df$sample)))+
##   theme(legend.position = "top")+ggtitle("ChIP-seq SET cross corr")
## print(p);dev.off()
