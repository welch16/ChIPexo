
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(spp)
library(parallel)
library(knitr)
library(reshape2)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/chip.seq.set.RData")
load("data/sample.summary.RData")
source("R/PBC_functions.R")

figsdir = "inst/figs/pbc"

exo_pbc = do.call(rbind,lapply(exo,PBC))
pet_pbc_bad= do.call(rbind,lapply(pet,PBC))
pet_pbc = do.call(rbind,lapply(pet,PBC,TRUE))
set_pbc = do.call(rbind,lapply(set,PBC))

save(list = c("exo_pbc","pet_pbc_bad","pet_pbc","set_pbc"),file = "data/PCR_bottleneck_coeff.RData")

exo_df = pbc_to_df(exo_pbc)
pet_df = pbc_to_df(pet_pbc)
pet_df_bad= pbc_to_df(pet_pbc_bad)
set_df = pbc_to_df(set_pbc)

exo_df$seq = "exo"
pet_df$seq = "pet"
pet_df_bad$seq = "pet-bad"
set_df$seq = "set"

df = rbind(exo_df,pet_df,pet_df_bad,set_df)
df$set = factor(df$set)
df$seq = factor(df$seq)

pdf(file = file.path(figsdir,"PBC_boxplot.pdf"))
p = ggplot(subset(df,subset =type== "all" & seq!="pet-bad" ),aes(seq,pbc,fill = seq))+geom_boxplot()+
  scale_y_continuous(limits = c(0,1))
print(p)
dev.off()

pdf(file = file.path(figsdir,"PBC_boxplot_strand_by_set.pdf"))
p = ggplot(subset(df,subset = seq == "exo"),aes(type,pbc,fill = type))+geom_boxplot()+ggtitle("ChIP-exo PBC")+
  scale_y_continuous(limits = c(0,1))
print(p)
p = ggplot(subset(df,subset = seq == "pet"),aes(type,pbc,fill = type))+geom_boxplot()+ggtitle("ChIP-seq PET PBC")+
  scale_y_continuous(limits = c(0,1))
print(p)  
p = ggplot(subset(df,subset = seq == "set"),aes(type,pbc,fill = type))+geom_boxplot()+ggtitle("ChIP-seq SET PBC")+
  scale_y_continuous(limits = c(0,1))
print(p)  
dev.off()
