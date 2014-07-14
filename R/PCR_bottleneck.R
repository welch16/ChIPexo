

rm(list = ls())
library(knitr)
library(parallel)
library(GenomicAlignments)
library(spp)


dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder



PBC <- function(file,cap = Inf)
{
  require(GenomicAlignments)
  message(file)
  rr = readGAlignmentsFromBam(file,param = NULL)
  ss1 = subset(rr,subset =strand(rr)=="+")
  ss2 = subset(rr,subset = strand(rr)=="-")
  tab1 = table(start(ss1))
  tab2 = table(end(ss2))
  tab = table(c(start(ss1),end(ss2)))
  PBC1 = round(sum(tab1 == 1)/sum(tab1 >= 1 & tab1 <= cap),4)
  PBC2 = round(sum(tab2 == 1)/sum(tab2 >= 1 & tab2 <= cap),4)
  chip.data = read.bam.tags(file)
  table.chip.data = lapply(chip.data$tags,table)
  nUniq = sum(sapply(table.chip.data,function(x) sum(x==1)))
  nTotal =  sum(sapply(table.chip.data,length))
  PBC = round(nUniq/nTotal ,4)
  message("PBC +:",PBC1)
  message("PBC -:",PBC2)
  message("PBC:",PBC)
  return(c(PBC_plus=PBC1,PBC_minus=PBC2,PBC=PBC))  
}
 

ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)

cap = Inf

suppressWarnings(SET <- mclapply(ff[[3]],FUN = PBC,cap=cap,mc.cores =8))
names(SET) = files[[3]]
suppressWarnings(PET <- mclapply(ff[[2]],FUN = PBC,cap=cap,mc.cores = 8))
names(PET) = files[[2]]
suppressWarnings(EXO <- mclapply(ff[[1]],FUN = PBC,cap=cap,mc.cores =8))
names(EXO) = files[[1]]

save(list = c("SET","PET","EXO"),file = "../RData/PBC.RData")





