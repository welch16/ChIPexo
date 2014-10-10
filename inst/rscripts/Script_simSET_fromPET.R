

source("R/density_functions.R")

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(data.table)

load("data/chip.seq.pet.RData")

set.seed('1234')
sample_set_from_pet <- function(pet_gr)
{  
  idd = sort(elementMetadata(pet_gr)[["qname"]],index.return=TRUE)$ix
  pet_gr = pet_gr[idd]
  nreads = max(elementMetadata(pet_gr)[["qname"]])
  bernoulli_sample = ifelse(runif(nreads)<=.5,1,2)
  idx = seq(0,nreads-1)*2 + bernoulli_sample
  return(pet_gr[idx])
}

set = mclapply(pet,sample_set_from_pet,mc.cores=8)
save(list = "set",file = "data/chip.seq.set.from.pet.RData")

