
rm(list = ls())
library(GenomicAlignments)
library(reshape2)
library(ggplot2)
library(parallel)

datadir = "data"
readfiles = c("Sig70_0min_sep_reads.RData",
  "Sig70_20min_sep_reads.RData",
  "BetaPrimeF_0min_sep_reads.RData",
  "BetaPrimeF_20min_sep_reads.RData")

files = c("outputSig70_0min_sep_reads.RData",
  "outputSig70_20min_sep_reads.RData",
  "outputBetaPrimeF_0min_sep_reads.RData",
  "outputBetaPrimeF_20min_sep_reads.RData")

load_out <- function(file)
{
  load(file)
  l1 = list(regions,exo1,exo2,pet1,pet2,set1,set2)
  names(l1) = c("regions","exo1","exo2","pet1","pet2","set1","set2")
  return(l1)
}

create_df <- function(alldata,ip,rif)
{
  regions = alldata[[1]]
  n = length(regions)
  df = elementMetadata(regions)
  df$width = width(regions)
  df$ip = ip
  df$rif = rif  
  return(df)
}


positions <- function(i,readsList)
{
  reads = readsList[[i]]
  if(is.null(reads)){
    out = 0
  }else{
    fwd = subset(reads,strand(reads) =="+")
    bwd = subset(reads,strand(reads) =="-")
    out = length(unique(start(fwd))) + length(unique(end(bwd)))    
  }
  return(out)
}


ip = rep(c("Sig70","BetaPrimeF"),each=2)
rif = rep(c("0min","20min"),2)


alldata = mclapply(file.path(datadir,files),load_out,mc.cores =4)
baseDf = mcmapply(create_df,alldata,ip,rif,SIMPLIFY=FALSE,mc.cores=4)


posList = list()

for(j in 1:4){
  load(file = file.path(datadir,readfiles[j]))
  n = nrow(baseDf[[j]])
  idxs = as.character(1:n)
  pos1 =  do.call(c,mclapply(idxs,positions,exo_sep_reads1,mc.cores=8))
  pos2 =  do.call(c,mclapply(idxs,positions,exo_sep_reads2,mc.cores=8))
  out = list(pos1,pos2)
  posList[[j]] = out
  
}


save(list = "posList",file = file.path(datadir,"position_numer.RData"))




