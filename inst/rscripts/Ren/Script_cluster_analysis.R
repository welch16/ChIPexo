# Read cluster analysis for Ren's data

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
#library(RColorBrewer)
#library(scales)

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

ourRegions = mcmapply(create_regions,regions_fwd,regions_bwd,MoreArgs = list(minW=0),SIMPLIFY=FALSE,mc.cores =6)
chr = paste0("chr",c(1:22,"X","Y"))

separate_regions <- function(regions,chr)
{
  regions = mclapply(chr,function(x,regions){
    u = subset(regions,seqnames == x)
    names(u) = NULL
    return(u)},regions,mc.cores =8)
  regions = GRangesList(regions)
  return(regions)
}

ourRegions = lapply(ourRegions,separate_regions,chr)
ourReads = lapply(reads,separate_regions,chr)


## dr = "../RenData/BAMfiles"


## # AY552.R1.
## param = ScanBamParam(which = ourRegions[[1]])
## reads = readGAlignmentsFromBam(file.path(dr,names(ourRegions)[1]),param = param,use.names = FALSE)

## save(list = "ourRegions",file = file.path(datadir,"Ren_regions.RData"))

## load(file =file.path(datadir,paste0("Ren_",names(ourRegions)[1],"_sepReads.RData")))


separate_reads<- function(reads,regions)
{
  overlaps = findOverlaps(regions,reads)
  sep_index = split(overlaps,factor(queryHits(overlaps)))
  separated_reads <- mclapply(sep_index,
    function(x,reads)reads[subjectHits(x)],reads,mc.cores = 8)
  return(separated_reads)
}



ourReads_sep = list()
j=5
  print(j)
  ourReads_sep[[j]] = mapply(separate_reads,ourReads[[j]],ourRegions[[j]],SIMPLIFY = FALSE)
  sep_reads = ourReads_sep[[j]]
  reg = ourRegions[[j]]
  save(list = c("sep_reads","reg"),file =file.path(datadir,paste0("Ren_",names(ourReads)[j],"_sepReads.RData")))


## names(ourReads_sep) = names(ourReads)
## save(file = file.path(datadir,"Ren_separated_reads.RData"),list = c("ourReads_sep","ourRegions"))


