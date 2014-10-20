
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(wavelets)
library(RColorBrewer)
library(scales)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/chip.seq.set.from.pet.RData")
load("data/sample.summary.RData")

source("R/depth_functions.R")
source("R/density_functions.R")
source("R/structure_functions.R")

# Parameters
mc.cores = 8
distance = 500

# Set conditions
conditions = rep("",4)
conditions[1] = resume.samples(ip="Sig70",rif="0 min")
conditions[2] = resume.samples(ip="BetaPrimeFlag",rif="0 min")
conditions[3] = resume.samples(ip="Sig70",rif="20 min")
conditions[4] = resume.samples(ip="BetaPrimeFlag",rif="20 min")

edsn = lapply(conditions,function(x,sample.info)as.character(subset(sample.info,
  subset = eval(parse(text=x)))$edsn),sample.info)
names(edsn) =conditions

exo.sets = lapply(edsn,function(x)exo[do.call(c,lapply(x,function(y)grep(y,names(exo))))])
pet.sets = lapply(edsn,function(x)pet[do.call(c,lapply(x,function(y)grep(y,names(pet))))])
set.sets = lapply(edsn,function(x)pet[do.call(c,lapply(x,function(y)grep(y,names(set))))])


exo.sets = lapply(exo.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)
pet.sets = lapply(pet.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)
set.sets = lapply(set.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)

indexdir = "inst/extdata"
indexfiles = c("Sig70_rif0_peak_plots_withIndicator.csv",
               "BetaPrimeFlag_rif0_peak_plots_withIndicator.csv",
               "Sig70_rif20_peak_plots_withIndicator.csv",
               "BetaPrimeFlag_rif20_peak_plots_withIndicator.csv")

extractFromCSV <- function(indexfile)
{
  table = read.table(indexfile,header =TRUE,sep =",") )
  ranges = IRanges(start = table$start, end = table$end)
  out = list(ranges,table$isPeak)
  return(out)
}

genomicRanges_baseTable <- function(indexfile,seqlength = seqlengths(exo.sets[[1]][[1]]))
{
  rangesList = extractFromCSV(indexfile)
  gr = GRanges(seqnames = names(seqlength),ranges =rangesList[[1]],strand = "*")
  gr$isPeak = rangesList[[2]]
  return(gr)
}


load("data/indexes_regions_withSet.RData")

## distance =500
## indexes = mapply(FUN = filterSets,exo.sets,pet.sets,MoreArgs = list(distance),SIMPLIFY = FALSE)
## names(indexes) = conditions


summary_stats_gr = lapply(file.path(indexdir,indexfiles),FUN = genomicRanges_baseTable)#,mc.cores=mc.cores)
names(summary_stats_gr) = conditions

exo.depths = lapply(exo.sets,function(x)lapply(x,length))
pet.depths = lapply(pet.sets,function(x)lapply(x,function(y)length(y)/2))
set.depths = lapply(set.sets,function(x)lapply(x,length))



regionReads <- function(summary_table,index_cond,seqSet,rep)
{  
  colname = paste0(seqSet,rep)
  elementMetadata(summary_table)[[colname]] = sapply(index_cond[[seqSet]][[rep]],length)
  return(summary_table)
}

summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("exo",1))
summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("exo",2))
summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("pet",1))
summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("pet",2))
summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("set",1))
summary_stats_gr = mapply(regionReads,summary_stats_gr,indexes,MoreArgs = list("set",2))


regionReadsScaled <- function(summary_table,depths,seqSet,rep)
{  
  colname = paste0(seqSet,rep)
  newColname = paste0(colname,"scaled")
  elementMetadata(summary_table)[[newColname]] = round(elementMetadata(summary_table)[[colname]]/depths[[rep]],6)
  return(summary_table) 
}

summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,exo.depths,MoreArgs = list("exo",1))
summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,exo.depths,MoreArgs = list("exo",2))
summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,pet.depths,MoreArgs = list("pet",1))
summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,pet.depths,MoreArgs = list("pet",2))
summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,pet.depths,MoreArgs = list("set",1))
summary_stats_gr = mapply(regionReadsScaled,summary_stats_gr,pet.depths,MoreArgs = list("set",2))



exo.reads = list()
pet.reads = list()
set.reads = list()
for(i in 1:4)
{
  exo.reads[[i]] = mapply(FUN = getReads,indexes[[i]]$exo,exo.sets[[i]],MoreArgs = list(mc.cores),SIMPLIFY = FALSE)
  pet.reads[[i]] = mapply(FUN = getReads,indexes[[i]]$pet,pet.sets[[i]],MoreArgs = list(mc.cores),SIMPLIFY = FALSE)
  set.reads[[i]] = mapply(FUN = getReads,indexes[[i]]$set,set.sets[[i]],MoreArgs = list(mc.cores),SIMPLIFY = FALSE)
}
names(exo.reads) = conditions
names(pet.reads) = conditions
names(set.reads) = conditions


strand_depth <- function(summary_table,reads,seqSet,rep,mc.cores)
{
  fwd_length = unlist(mclapply(reads[[rep]],function(x,str)length(subset(x,strand == str)),"+",mc.cores =mc.cores))
  bwd_length = unlist(mclapply(reads[[rep]],function(x,str)length(subset(x,strand == str)),"-",mc.cores =mc.cores))
  colname = paste0(seqSet,rep)
  elementMetadata(summary_table)[[paste0(colname,"_fwdDepth")]] = fwd_length
  elementMetadata(summary_table)[[paste0(colname,"_bwdDepth")]] = bwd_length
  return(summary_table)  
}

summary_stats_gr = mapply(strand_depth,summary_stats_gr,exo.reads,MoreArgs = list("exo",1,mc.cores))
summary_stats_gr = mapply(strand_depth,summary_stats_gr,exo.reads,MoreArgs = list("exo",2,mc.cores))
summary_stats_gr = mapply(strand_depth,summary_stats_gr,pet.reads,MoreArgs = list("pet",1,mc.cores))
summary_stats_gr = mapply(strand_depth,summary_stats_gr,pet.reads,MoreArgs = list("pet",2,mc.cores))
summary_stats_gr = mapply(strand_depth,summary_stats_gr,pet.reads,MoreArgs = list("set",1,mc.cores))
summary_stats_gr = mapply(strand_depth,summary_stats_gr,pet.reads,MoreArgs = list("set",2,mc.cores))



fwd_strand_ratio <- function(summary_table,reads,seqSet,rep)
{
  colname1 = paste0(seqSet,rep,"_fwdDepth")
  colname2 = paste0(seqSet,rep,"_bwdDepth")
  f = elementMetadata(summary_table)[[colname1]]
  r = elementMetadata(summary_table)[[colname2]]
  elementMetadata(summary_table)[[paste0(seqSet,rep,"_fwdStrandRatio")]] = round(ifelse(f+r==0,0,f/(f+r)),4)
  return(summary_table)
}

summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,exo.reads,MoreArgs = list("exo",1))
summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,exo.reads,MoreArgs = list("exo",2))
summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,pet.reads,MoreArgs = list("pet",1))
summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,pet.reads,MoreArgs = list("pet",2))
summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,pet.reads,MoreArgs = list("set",1))
summary_stats_gr = mapply(fwd_strand_ratio,summary_stats_gr,pet.reads,MoreArgs = list("set",2))


summit_from_cover <- function(cover,lb,ub)
{
  out = NA
  if(nrun(cover)>1){
    out = which.max(coverToVec(lb,ub,cover))
  }
  return(out)
}

strand_summit_position <- function(summary_table,reads,seqSet,rep,mc.cores)
{
  cover_F = mclapply(reads[[rep]],function(x,str)cover_by_strand(x,str),"+",mc.cores =mc.cores)
  cover_R = mclapply(reads[[rep]],function(x,str)cover_by_strand(x,str),"-",mc.cores =mc.cores)
  summit_F = unlist(mcmapply(summit_from_cover,cover_F,start(summary_table),end(summary_table),
    mc.cores =mc.cores,SIMPLIFY=FALSE))
  summit_R = unlist(mcmapply(summit_from_cover,cover_R,start(summary_table),end(summary_table),
    mc.cores =mc.cores,SIMPLIFY=FALSE))
  diff = summit_F - summit_R
  elementMetadata(summary_table)[[paste0(seqSet,rep,"_fwdSummit_pos")]] = summit_F
  elementMetadata(summary_table)[[paste0(seqSet,rep,"_bwdSummit_pos")]] = summit_R
  elementMetadata(summary_table)[[paste0(seqSet,rep,"_summit_pos_diff")]] = diff
  return(summary_table)
}

summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,exo.reads,MoreArgs = list("exo",1,mc.cores))
summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,exo.reads,MoreArgs = list("exo",2,mc.cores))
summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,pet.reads,MoreArgs = list("pet",1,mc.cores))
summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,pet.reads,MoreArgs = list("pet",2,mc.cores))
summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,pet.reads,MoreArgs = list("set",1,mc.cores))
summary_stats_gr = mapply(strand_summit_position,summary_stats_gr,pet.reads,MoreArgs = list("set",2,mc.cores))



save(list = "summary_stats_gr",file = "data/summary_stats_labeled_regions_withSet.RData")


