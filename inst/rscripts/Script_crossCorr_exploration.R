
rm(list = ls())
library(parallel)
library(spp)
library(ChIPQC)
library(reshape2)
library(ggplot2)
source("R/cross_corr_functions.R")

# This script only converts the data from bam file format to RData

dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")
range = c(1,300)
bin = 1

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder
ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)
names(ff) = folder

# ma <- function(x,n)filter(x,rep(1/n,n),sides = 2)

mc.cores = 8
exo_chipqc = read_chipqc(ff[[1]],mc.cores)
pet_chipqc = read_chipqc(ff[[2]],mc.cores)
set_chipqc = read_chipqc(ff[[3]],mc.cores)

exo_spp = read_spp(ff[[1]],mc.cores)
pet_spp = read_spp(ff[[2]],mc.cores)
set_spp = read_spp(ff[[3]],mc.cores)

exo_crossCorr_spp = cross_corr_spp(exo_spp,range,bin)
pet_crossCorr_spp = cross_corr_spp(pet_spp,range,bin)
set_crossCorr_spp = cross_corr_spp(set_spp,range,bin)

shift = seq(range[1],range[2],by=bin)


exo_df_list = cross_corr_df(exo_chipqc,exo_crossCorr_spp,files[[1]],shift,mc.cores)
pet_df_list = cross_corr_df(pet_chipqc,pet_crossCorr_spp,files[[2]],shift,mc.cores)
set_df_list = cross_corr_df(set_chipqc,set_crossCorr_spp,files[[3]],shift,mc.cores)


# load("data/cross_corr.RData")

nrow = length(exo_df_list)
pdf(file = "inst/figs/cross_corr/ChIPExo_Cross_Corr.pdf",height = 2*nrow)
p = plot_template(exo_df_list,"ChIP-Exo Cross Corr",nrow)
print(p);dev.off()

nrow = length(pet_df_list)
pdf(file = "inst/figs/cross_corr/ChIPSeq_PET_Cross_Corr.pdf",height = 2*nrow)
p = plot_template(pet_df_list,"ChIP-Seq PET Cross Corr",nrow)
print(p);dev.off()

nrow = length(set_df_list)
pdf(file = "inst/figs/cross_corr/ChIPSeq_SET_Cross_Corr.pdf",height = 2*nrow)
p = plot_template(set_df_list,"ChIP-Seq SET Cross Corr",nrow)
print(p);dev.off()

save(list = c("exo_df_list","pet_df_list","set_df_list"),file = "data/cross_corr.RData")
