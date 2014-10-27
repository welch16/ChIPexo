#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_slice_summary_stat.R - Adds a list of summary statistics to a set of regions

  Arguments:

    -- regionsFile

      The RData file to be updated

   -- help

      Show the help file

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

if(length(args)!=1){
  stop("--regionsFile must be specified")
}


regionsFile = args[1]

datadir = "data"
load(file.path(datadir,"regions_slice.RData"))



if(regionsFile == "Sig70_0min_sep_reads.RData"){
  regions = ourRegions[[1]]
}else if(regionsFile == "Sig70_20min_sep_reads.RData"){
  regions = ourRegions[[3]]
}else if(regionsFile == "BetaPrimeF_0min_sep_reads.RData"){
  regions = ourRegions[[2]]
}else if(regionsFile == "BetaPrimeF_20min_sep_reads.RData"){
  regions = ourRegions[[4]]
}else{
  message("Wrong....\n")
  q()}

library(GenomicAlignments)



load(file.path(datadir,regionsFile))

diff_pos_summit <- function(reads)
{
  fwd_cover = suppressWarnings(coverage(subset(reads,strand(reads)=="+"))[[1]])
  bwd_cover = suppressWarnings(coverage(subset(reads,strand(reads)=="-"))[[1]])
  max_fwd = NA
  fwd_summit_pos = NA
  if(nrun(fwd_cover) > 1){
    max_fwd = max(fwd_cover)
    fwd_summit_pos = head(which(fwd_cover == max_fwd),n=1)
  }
  max_bwd = NA
  bwd_summit_pos = NA
  if(nrun(bwd_cover) > 1){
    max_bwd = max(bwd_cover)
    bwd_summit_pos = tail(which(bwd_cover == max_bwd),n=1)
  }  
  diff = bwd_summit_pos - fwd_summit_pos    
  return(c(diff=diff,max_fwd=max_fwd,max_bwd=max_bwd,fwd_pos = fwd_summit_pos,bwd_pos = bwd_summit_pos))
}

fwd_strand_ratio <- function(reads)
{
  f = suppressWarnings(length(subset(reads,strand(reads)=="+")))
  r = suppressWarnings(length(subset(reads,strand(reads)=="-")))
  if( f+r == 0){
    prob = NA
  }else{
    prob = f / (f+r)
  }
  return(c(prob=prob,f=f,r=r))
}

library(parallel)

summary_statistics <- function(reads,regions,mc)
{
  fwd_strand = t(do.call(cbind,mclapply(1:length(regions),function(i,reads)fwd_strand_ratio(reads[[i]]),reads,mc.cores = mc)))
  dif_pos = t(do.call(cbind,mclapply(1:length(regions),function(i)diff_pos_summit(reads[[i]]),mc.cores = mc)))
  df = DataFrame(cbind(fwd_strand,dif_pos))
  return(df)
}
  
# Exo1
exo1 = summary_statistics(exo_sep_reads1,regions,8)
# Exo2
exo2 = summary_statistics(exo_sep_reads2,regions,8)
# Pet1
pet1 = summary_statistics(pet_sep_reads1,regions,8)
# Pet2
pet2 = summary_statistics(pet_sep_reads2,regions,8)
# Set1
set1 = summary_statistics(set_sep_reads1,regions,8)
# Set2
set2 = summary_statistics(set_sep_reads2,regions,8)

save(file = file.path(datadir,regionsFile),list = c("regions","exo1","exo2","pet1","pet2","set1","set2"))
  
q()
