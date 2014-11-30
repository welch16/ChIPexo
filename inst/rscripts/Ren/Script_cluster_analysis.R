# Read cluster analysis for Ren's data

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(pryr)
library(plyr)

datadir = "data"
load(file = file.path(datadir,"Ren_reads.RData"))

readsR1 = list(reads1[[1]],reads2[[1]],reads3[[1]])
names(readsR1) = paste0("AY55",2:4)

## readsR2 = list(reads1[[2]],reads2[[2]],reads3[[2]])

reads_fwd = mclapply(readsR1,function(x)subset(x,strand(x)=="+"),mc.cores= 3)
reads_bwd = mclapply(readsR1,function(x)subset(x,strand(x)=="-"),mc.cores= 3)

cover_fwd = mclapply(reads_fwd,coverage,mc.cores = 3)
cover_bwd = mclapply(reads_bwd,coverage,mc.cores = 3)

slice_all <- function(cover,q)
{
  lw = quantile(runValue(cover),q)
  slice_list= mcmapply(slice,cover,lw,MoreArgs = list(rangesOnly=TRUE),SIMPLIFY=FALSE,mc.cores=8)
  slice_list = GRangesList(mcmapply(function(x,y)GRanges(seqnames = y,ranges = x,strand = "*"),slice_list,names(slice_list),
    SIMPLIFY= FALSE,mc.cores = 8))
  return(unlist(slice_list))
}

q = 0.75

regions_fwd = lapply(cover_fwd,slice_all,q=q)
regions_bwd = lapply(cover_bwd,slice_all,q=q)

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
  names(out) = NULL
  return(out)
}

ourRegions = mcmapply(create_regions,regions_fwd,regions_bwd,MoreArgs = list(minW=0),SIMPLIFY=FALSE,mc.cores =3)
readLength = 40

# Filter the regions with width less than readLength
ourRegions = lapply(ourRegions,function(x)x[width(x)>readLength])
chr = paste0("chr",c(1:22,"X","Y"))

separate_regions <- function(regions,chr)
{
  regions = mclapply(chr,function(x,regions){
    u = subset(regions,seqnames == x)
    u = u[order(u)]
    return(u)},regions,mc.cores =8)
  regions = GRangesList(regions)
  names(regions) = chr
  return(regions)
}

regions = list()
regions[["both"]] = lapply(ourRegions,function(x,lb)
  subset(x,label == lb),"both")
regions[["fwd"]] = lapply(ourRegions,function(x,lb)
  subset(x,label == lb),"fwd")
regions[["bwd"]] = lapply(ourRegions,function(x,lb)
  subset(x,label == lb),"bwd")

regions[["both"]] = lapply(regions[["both"]],separate_regions,chr)
regions[["fwd"]] = lapply(regions[["fwd"]],separate_regions,chr)
regions[["bwd"]] = lapply(regions[["bwd"]],separate_regions,chr)

reads = list()
reads[["fwd"]] = lapply(reads_fwd,separate_regions,chr)
reads[["bwd"]] = lapply(reads_bwd,separate_regions,chr)

dr = "/p/keles/ChIPexo/volume3/ChIPexo/data/Ren/separated"

save(list = "reads",file = file.path(dr,"separated_reads.RData"))
save(list = "regions",file = file.path(dr,"separated_regions.RData"))

overlaps = list()
sets = names(reads[[1]])
for(set in sets){
  overlaps[[set]] = list(
    both = list(
      fwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["both"]][[set]],
        reads[["fwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12),
      bwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["both"]][[set]],
        reads[["bwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12)
    ),
    fwd = list(
      fwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["fwd"]][[set]],
        reads[["fwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12),
      bwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["fwd"]][[set]],
        reads[["bwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12)
    ),
    bwd = list(
      fwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["bwd"]][[set]],
        reads[["fwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12),
      bwd = mcmapply(function(x,y)findOverlaps(x,y),
        regions[["bwd"]][[set]],
        reads[["bwd"]][[set]],SIMPLIFY=FALSE,mc.cores= 12)
    )            
  )
}

save(list = "overlaps",file = file.path(dr,"separated_overlaps.RData"))




## separate_reads<- function(reads,regions)
## {
##   # can make it smarter, dont need to look fwd or bwd
##   overlaps = findOverlaps(regions,reads)
##   sep_index = split(overlaps,factor(queryHits(overlaps)))
##   browser()
##   separated_reads <- mclapply(sep_index,
##     function(x,reads)reads[subjectHits(x)],reads,mc.cores = 8)
##   return(separated_reads)
## }





## separate_reads(reads_fwd[[1]][[23]],regions_both[[1]][[23]])


## u1 = uu[[24]]
## v1 = vv[[24]]
## our1 = ourRegions[[1]][[24]]
## our1 = our1[order(our1)]

## pp = separate_reads(u1,our1)



## ourReads_sep = list()
## j=5
##   print(j)
##   ourReads_sep[[j]] = mapply(separate_reads,ourReads[[j]],ourRegions[[j]],SIMPLIFY = FALSE)
##   sep_reads = ourReads_sep[[j]]
##   reg = ourRegions[[j]]
##   save(list = c("sep_reads","reg"),file =file.path(datadir,paste0("Ren_",names(ourReads)[j],"_sepReads.RData")))


## names(ourReads_sep) = names(ourReads)
## save(file = file.path(datadir,"Ren_separated_reads.RData"),list = c("ourReads_sep","ourRegions"))


