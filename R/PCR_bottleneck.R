

rm(list = ls())


library(knitr)
library(parallel)
library(GenomicAlignments)

dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder



PBC <- function(file)
{
  require(GenomicAlignments)
  message(file)
  rr = readGAlignmentsFromBam(file,param = NULL)
  ss1 = subset(rr,subset =strand(rr)=="+")
  ss2 = subset(rr,subset = strand(rr)=="-")
  N11 = length(unique(start(ss1)))
  N12 = length(unique(end(ss2)))
  N1 = N11+N12
  Nd1 = length(ss1)
  Nd2 = length(ss2)
  Nd = Nd1 + Nd2
  PBC1 = round(N11 / Nd1,4)
  message("PBC +:",PBC1)
  PBC2 = round(N12 / Nd2,4)
  message("PBC -:",PBC2)
  PBC = round(N1 /Nd,4)
  message("PBC:",PBC)
  return(c(PBC_plus=PBC1,PBC_minus=PBC2,PBC=PBC))  
}
 

ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)

suppressWarnings(SET <- lapply(ff[[3]],FUN = PBC))
names(SET) = files[[3]]
suppressWarnings(PET <- lapply(ff[[2]],FUN = PBC))
names(PET) = files[[2]]
suppressWarnings(EXO <- lapply(ff[[1]],FUN = PBC))
names(EXO) = files[[1]]

save(list = c("SET","PET","EXO"),file = "../RData/PBC.RData")
