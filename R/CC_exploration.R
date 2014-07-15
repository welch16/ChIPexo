
rm(list = ls())
library(parallel)
library(spp)
library(snow)
library(GenomicAlignments)
library(ChIPQC)


dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder
ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)
names(ff) = folder

ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)

nsc_matrix <- function(cc)
{
  # cc is a data frame with 2 columns: shift and cc_score  
  shift = cc[,1]
  cc_score = cc[,2]
  max_shift = length(shift)
  mat = matrix(NA,max_shift,max_shift)
  for(s in 1:max_shift){
    for(t in 1:s){
      mat[s,t] = cc_score[s] / cc_score[t]
    }
  }
  return(mat)
}

read_data <- function(file)
{
  qc.sample = ChIPQCsample(file,annotation = NULL)
  spp.sample = read.bam.tags(file)
  return(list(chipqc= qc.sample,spp=spp.sample))
}

exo.data = lapply(ff[[1]],FUN = read_data)
names(exo.data) = files[[1]]

pet.data = lapply(ff[[2]],FUN = read_data)
names(pet.data) = files[[2]]

set.data = lapply(ff[[3]],FUN = read_data)
names(set.data) = files[[3]]

save(list = c("exo.data","pet.data","set.data"),file = "bam_data_chipqc_spp.RData")

mc = 12

qc_metrics_cross_corr <- function(qc.sample)
{
  rf = QCmetrics(qc.sample)
  cc_chipqc= plotCC(qc.sample)$data
  return(list(metrics = rf,cross.correlation = cc_chipqc))
}

spp_binding_char <- function(spp.sample)
{
  return(get.binding.characteristics(chip.data,srange = c(1,300),bin = 1))
}

fragLen_and_cc <- function(seq.data)
{
  qc.metrics = qc_metrics_cross_corr(seq.data$chipqc)
  qc.readL = qc.metrics$metrics[5]
  qc.cross = qc.metrics$cross.correlation
  qc.fragL = qc.metrics$metrics[6]
  spp.metrics = spp_binding_char(seq.data$spp)
  spp.fragL = spp.metrics$peak$x
  spp.cross = spp.metrics$cross.correlation
  return(list(qc.readL=qc.readL,qc.fragL = qc.fragL,spp.fragL = spp.fragL,qc.cross = qc.cross,spp.cross = spp.cross))
}

exo.cc = mclapply(exo.data,FUN = fragLen_and_cc,mc.cores = mc)
names(exo.cc) = files[[1]]

pet.cc = mclapply(pet.data,FUN = fragLen_and_cc,mc.cores = mc)
names(pet.cc) = files[[2]]

set.cc = mclapply(set.data,FUN = fragLen_and_cc,mc.cores = mc)
names(set.cc) = files[[3]]

save(list = c("exo.cc","pet.cc","set.cc"),file = "chipqc_spp_fragL_cros_corr.RData")




