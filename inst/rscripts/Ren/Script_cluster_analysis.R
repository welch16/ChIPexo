# Read cluster analysis for Ren's data

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(wavelets)
library(RColorBrewer)
library(scales)

datadir = "data"
load(file = file.path(datadir,"Ren_reads.RData"))

reads = c(reads1,reads2,reads3)

reads_fwd = mclapply(reads,function(x)subset(x,strand(x)=="+"),mc.cores= 6)
reads_bwd = mclapply(reads,function(x)subset(x,strand(x)=="-"),mc.cores= 6)

cover_fwd = mclapply(reads_fwd,coverage,mc.cores = 6)
cover_bwd = mclapply(reads_bwd,coverage,mc.cores = 6)


slice_all <- function(cover,lw)
{
  slice_list= slice(cover,lower = lw,rangesOnly=TRUE)
  slice_list = GRangesList(mcmapply(function(x,y)GRanges(seqnames = y,ranges = x,strand = "*"),slice_list,names(slice_list),
    SIMPLIFY= FALSE,mc.cores = 8))
  return(unlist(slice_list))
}

regions_fwd = lapply(cover_fwd,slice_all,lw=5)
regions_bwd = lapply(cover_bwd,slice_all,lw=5)

create_regions <- function(fwd,bwd,minW =0)
{
  fwd_l = length(fwd)
  bwd_l = length(bwd)
  overlaps = findOverlaps(fwd,bwd)
  fwd_only = fwd[!(1:fwd_l %in% queryHits(overlaps))]
  bwd_only = bwd[!(1:bwd_l %in% subjectHits(overlaps))]
  both = union(fwd[queryHits(overlaps)],bwd[subjectHits(overlaps)])
  both$label = "both"
  fwd_only$label = "fwd"
  bwd_only$label = "bwd"
  out = c(both,fwd_only,bwd_only)
  out = out[width(out) > minW]
  return(out)
}

ourRegions = mcmapply(create_regions,regions_fwd,regions_bwd,MoreArgs = list(minW=50),SIMPLIFY=FALSE,mc.cores =6)


separated_reads<- function(reads,regions)
{
  overlaps = findOverlaps(regions,reads)
  n_regions = length(regions)
  separated_reads = mclapply(1:n_regions,function(i,overlaps,reads)reads[subjectHits(overlaps[queryHits(overlaps)==i])],overlaps,reads,mc.cores = 8)
  return(separated_reads)
}

sep_reads = mapply(separated_reads,reads,ourRegions,SIMPLIFY = FALSE)

save(file = file.path(datadir,"Ren_separated_reads.RData"),list = "sep_reads")
