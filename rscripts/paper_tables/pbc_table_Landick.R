
rm(list = ls())

source("R/base_edsn.R")

library(data.table)
library(GenomicAlignments)
library(ChIPUtils)

## load characteristics

what <- c("exo","pet","set")

char1 <- lapply(what,edsn_tab)
char2 <- lapply(what,edsn_tab_old)
char3 <- copy(char2[[1]])
char3[,edsn := "933"]

#######################################################################################

char1 <- lapply(char1,function(x)x[ip == "Sig70"])

#######################################################################################

### reads directories

chip_dirs <- list()
chip_dirs[["exo"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
chip_dirs[["pet"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
chip_dirs[["set"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

######################################################################################

## get files

read_files1 <- mapply(function(dirs,char,pet){
  cc <- split(char,char$edsn)
  out <- sapply(cc,function(x)filter_reads(dirs,x,pet))
  names(out) <- NULL
  return(out)},chip_dirs,char1,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

read_files2 <- mapply(filter_reads,
  chip_dirs,char2,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

read_files3 <- mapply(filter_reads,
  chip_dirs[1],list(char3),c(FALSE),SIMPLIFY = FALSE)

######################################################################################

reads2 <- mapply(create_reads,read_files2,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)
reads3 <- mapply(create_reads,read_files3,c(FALSE),SIMPLIFY = FALSE)

reads1 <- mapply(function(files,pet){
  out <- lapply(files,create_reads,pet)
  return(out)},read_files1,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)


pbc1 <- lapply(reads1,function(x)sapply(x,PBC))
pbc2 <- lapply(reads2,PBC)
pbc3 <- lapply(reads3,PBC)

char1 <- mapply(function(x,y)x[,pbc := y],char1,pbc1,SIMPLIFY = FALSE)
names(char1) <- what
char2 <- mapply(function(x,y)x[,pbc := y],char2,pbc2,SIMPLIFY = FALSE)
names(char2) <- what
char3[,pbc := pbc3[[1]]]

pbc <- list(char1,char2,char3)

save(pbc,file = "data/for_paper/EColi_experiments_PBC_tables.RData")
