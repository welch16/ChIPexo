

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/chip.seq.set.RData")
load("data/sample.summary.RData")

distance = 500
mc.cores = 8

exo = mclapply(exo,as.GRanges,mc.cores=mc.cores)
pet = mclapply(pet,as.GRanges,mc.cores=mc.cores)
set = mclapply(set,as.GRanges,mc.cores=mc.cores)


