

rm(list = ls())
source("R/PBC_functions.R")
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)


datadir = "data"
load(file = file.path(datadir,"Ren_reads.RData"))

reads = c(reads1,reads2,reads3)

pbc = lapply(reads,PBC)

save(list = "pbc",file = file.path(datadir,"PBC_RenData.RData"))

load(file.path(datadir,"PBC_RenData.RData"))
pbc = as.data.frame(t(do.call(cbind,pbc)))
pbc$adapter = rep(c("R1","R2"),3)
pbc$sample = rep(c("AY552","AY553","AY554"),each = 2)
figsdir = "inst/figs/pbc"
pdf(file.path(figsdir,"PBC_Ren.pdf"))
p1 = ggplot(pbc,aes(sample,all,colour = adapter))+geom_point(size = 3)
print(p1)
dev.off()
