
rm(list = ls())

library(GenomicAlignments)
library(parallel)
library(data.table)
library(pryr)

dr = "/p/keles/ChIPexo/volume3/ChIPexo/data/Ren/separated"

load(file = file.path(dr,"separated_reads_lw3.RData"))
# reads
load(file = file.path(dr,"separated_regions_lw3.RData"))
# regions
load(file = file.path(dr,"separated_overlaps_lw3.RData"))
# overlaps

# Want to for each data set obtain a set of summary statistics

sets = names(overlaps)
labels = names(overlaps[[1]])
strands = names(overlaps[[1]][[1]])
chr = names(overlaps[[1]][[1]][[1]])

extract_reads <- function(param,reads,regions,overlaps)
{
  ourReads = reads[[param[3]]][[param[1]]][[param[4]]]
  ourRegions = regions[[param[2]]][[param[1]]][[param[4]]]
  ourOverlaps = overlaps[[param[1]]][[param[2]]][[param[3]]][[param[4]]]
  ourReads = data.table(as.data.frame(ourReads))
  ourReads$ID = 1:nrow(ourReads)
  ourOverlaps = data.table(as.data.frame(ourOverlaps))
  setkey(ourOverlaps,queryHits)
  setkey(ourReads,ID)
  ourReads = ourReads[ourOverlaps[queryHits %in% 1:length(ourRegions),list(subjectHits)],]
  ourReads$region = ourOverlaps[,list(queryHits)]
  return(ourReads)
}


build_param_Df <- function(set,labels,strands,chr)
{
  paramDf = as.data.frame(matrix(nrow = length(labels)
    *length(strands)*length(chr),ncol=4))
  names(paramDf) = c("set","label","strand","chr")
  k=1
  for(lab in labels){
    for(stran in strands){
      for(ch in chr){
        paramDf[k,] = c(set,lab,stran,ch)
        k = k+1
      }
    }
  }
  paramDf = data.table(paramDf)
  paramDf$ID = 1:nrow(paramDf)
  return(paramDf)
}

paramDf1 = build_param_Df(sets[1],labels,strands,chr)
paramDf2 = build_param_Df(sets[2],labels,strands,chr)
paramDf3 = build_param_Df(sets[3],labels,strands,chr)

setkey(paramDf1,ID,label,strand,chr)
setkey(paramDf2,ID,label,strand,chr)
setkey(paramDf3,ID,label,strand,chr)


sep_reads1 = mclapply(1:nrow(paramDf1),function(i){
  param = as.character(paramDf1[ID==i,list(set,label,strand,chr)])
  separated_reads<- extract_reads(param,reads,regions,overlaps)
  return(separated_reads)},mc.cores = 16)

sep_reads2 = mclapply(1:nrow(paramDf2),function(i){
  param = as.character(paramDf2[ID==i,list(set,label,strand,chr)])
  separated_reads<- extract_reads(param,reads,regions,overlaps)
  return(separated_reads)},mc.cores = 16)

sep_reads3 = mclapply(1:nrow(paramDf3),function(i){
  param = as.character(paramDf3[ID==i,list(set,label,strand,chr)])
  separated_reads<- extract_reads(param,reads,regions,overlaps)
  return(separated_reads)},mc.cores = 16)


data.table2IRanges <- function(DT){
  ranges = IRanges(start = DT$start,
      end = DT$end)
  return(ranges)
}

data.table2GRanges <- function(DT){
  sn = as.character(DT$seqnames)
  gr = GRanges(seqnames = sn,
    ranges = data.table2IRanges(DT),
    strand = DT$strand)
  return(gr)
}


strand_statistics <- function(region_reads,param,idx)
{
  # depth
  depth = nrow(region_reads)
  # Nr. positions
  ranges= data.table2IRanges(region_reads)
  if(as.character(param[,list(strand)])=="fwd"){
    nrPos = length(unique(start(ranges)))
  }else{
    nrPos = length(unique(end(ranges)))
  }
  newStart = min(start(ranges))
  newEnd = max(end(ranges))
  # Summit position
  cover = coverage(ranges)
  mean = mean(window(cover,start=newStart,end =newEnd))
  sd = sd(window(cover,start=newStart,end =newEnd))
  snr = mean / sd
  if(nrun(cover) > 1){
    maxCover = max(cover)
    if(identical(param$strand,"bwd")){
      summitPos = max(which(cover == maxCover))
    }else{
      summitPos = min(which(cover == maxCover))
    }
  }else{
    maxCover = NA
    summitPos = NA
  }
  return(c(idx=idx,depth=depth,nrPos =nrPos,maxCover=maxCover,
           summitPos=summitPos,mean = mean,sd = sd,snr = snr,
           newstart  =newStart, newend = newEnd))
}

summaryStats1 = lapply(1:nrow(paramDf1),function(j){
  print(paste0(j,"/144"))
  reg = as.numeric(unique(sep_reads1[[j]][,list(region)])[[1]])
  stats = mclapply(reg,
    function(i)strand_statistics(sep_reads1[[j]][region==i,],paramDf1[ID==i,],i),mc.cores=24)
  return(stats)
})

summaryStats2 = lapply(1:nrow(paramDf2),function(j){
  print(paste0(j,"/144"))
  reg = as.numeric(unique(sep_reads2[[j]][,list(region)])[[1]])
  stats = mclapply(reg,
    function(i)strand_statistics(sep_reads2[[j]][region==i,],paramDf2[ID==i,],i),mc.cores=24)
  return(stats)
})


summaryStats3 = lapply(1:nrow(paramDf3),function(j){
  print(paste0(j,"/144"))  
  reg = as.numeric(unique(sep_reads3[[j]][,list(region)])[[1]])
  stats = mclapply(reg,
    function(i)strand_statistics(sep_reads3[[j]][region==i,],paramDf3[ID==i,],i),mc.cores=24)
  return(stats)
})

save(list = "summaryStats1",file = file.path(dr,"AY552_summary_lw3.RData"))
save(list = "summaryStats2",file = file.path(dr,"AY553_summary_lw3.RData"))
save(list = "summaryStats3",file = file.path(dr,"AY554_summary_lw3.RData"))

load(file = file.path(dr,"AY552_summary_lw3.RData"))
load(file = file.path(dr,"AY553_summary_lw3.RData"))
load(file = file.path(dr,"AY554_summary_lw3.RData"))

merge_case <- function(ourset,ourlabel,ourchr,summaryStats,paramDf,regions)
{
  ourRegions = regions[[ourlabel]][[ourset]][[ourchr]]
  ourParam = paramDf[label == ourlabel & set == ourset & chr == ourchr]
  ids = ourParam[,list(ID)][[1]]
  names(ids) = ourParam[,list(strand)][[1]]
  ourstats = summaryStats[ids]
  names(ourstats) = names(ids) # "fwd" & "bwd"
  if(nrow(ourstats[["bwd"]])==0){ # quick fix, need to fix this bug
    ourstats[["bwd"]] = ourstats[["bwd"]][idx==1,]
    ourstats[["bwd"]]$depth = 0
    ourstats[["bwd"]]$nrPos = 0
    ourstats[["bwd"]]$maxCover = NA
    ourstats[["bwd"]]$summitPos = NA
    ourstats[["bwd"]]$mean = NA
    ourstats[["bwd"]]$sd = NA
    ourstats[["bwd"]]$snr = NA
  }
  setkey(ourstats[["bwd"]],idx)
  if(nrow(ourstats[["fwd"]])==0){ # quick fix, need to fix this bug
    ourstats[["fwd"]] = ourstats[["bwd"]][idx==1,]
    ourstats[["fwd"]]$depth = 0
    ourstats[["fwd"]]$nrPos = 0
    ourstats[["fwd"]]$maxCover = NA
    ourstats[["fwd"]]$summitPos = NA
    ourstats[["fwd"]]$mean = NA
    ourstats[["fwd"]]$sd = NA
    ourstats[["fwd"]]$snr = NA    
  }
  setkey(ourstats[["fwd"]],idx)
  extract_summary <- function(sumstats){
    out = rep(NA,6)
    if(nrow(sumstats) > 0){
      out[1] = sumstats[["depth"]]
      out[2] = sumstats[["nrPos"]]
      out[3] = sumstats[["maxCover"]]      
      out[4] = sumstats[["summitPos"]]
      out[5] = sumstats[["newstart"]]
      out[6] = sumstats[["newend"]]
      out[7] = sumstats[["mean"]]
      out[8] = sumstats[["sd"]]
      out[9] = sumstats[["snr"]]
    }
    out[1:2] = ifelse(is.na(out[1:2]),0,out[1:2])
    return(out)    
  }
  ourStats = data.table(do.call(rbind,mclapply(1:length(ourRegions),function(i,ourstats){
    fwd = extract_summary(ourstats[["fwd"]][idx==i])
    bwd = extract_summary(ourstats[["bwd"]][idx==i])
    depth = fwd[1] + bwd[1]
    nr_positions = fwd[2] + bwd[2]
    prob = fwd[1] / depth
    diff = bwd[4] - fwd[4]
    newstart = min(fwd[5],bwd[5],na.rm=TRUE)
    newend = max(fwd[6],bwd[6],na.rm =TRUE)
    out = c(fwd[1:4],bwd[1:4],depth,nr_positions,prob,diff,newstart,newend)
    names(out) = c("f","f_nrpos","f_maxcover","f_position",
           "r","r_nrpos","r_maxcover","r_position",
           "depth","nr_pos","prob","diff","newstart","newend")
    return(out)    
  },ourstats,mc.cores =24)))
  return(ourStats)  
}

merge_stats <- function(paramDf,summaryStats,regions)
{  
  summaryStats = mclapply(summaryStats,function(x)do.call(rbind,x),mc.cores=24)
  summaryStats = mclapply(summaryStats,function(x)data.table(x),mc.cores=24)
  setkey(paramDf,strand)
  subParamDf = paramDf[strand == "fwd"]
  subParamDf$idx = 1:nrow(subParamDf)
  setkey(subParamDf,idx)  
  merged_stats = lapply(1:nrow(subParamDf),function(j){
    print(paste0(j,"/72"))
    ourparam = subParamDf[idx == j]
    merged_stats = merge_case(ourparam[["set"]],ourparam[["label"]],
      ourparam[["chr"]],summaryStats,paramDf,regions)
    return(merged_stats)})
  return(merged_stats)
}

merged_stats1 = merge_stats(paramDf1,summaryStats1,regions)
merged_stats2 = merge_stats(paramDf2,summaryStats2,regions)
merged_stats3 = merge_stats(paramDf3,summaryStats3,regions)

save(list = "merged_stats1",file = file.path(dr,"AY552_summary_merged_lw3.RData"))
save(list = "merged_stats2",file = file.path(dr,"AY553_summary_merged_lw3.RData"))
save(list = "merged_stats3",file = file.path(dr,"AY554_summary_merged_lw3.RData"))

## load(file = file.path(dr,"AY552_summary_merged.RData"))
## load(file = file.path(dr,"AY553_summary_merged.RData"))
## load(file = file.path(dr,"AY554_summary_merged.RData"))

merge_all <- function(merged_stats,paramDf,regions)
{
  labels = paramDf[,list(label)][[1]]
  chrs = paramDf[,list(chr)][[1]]
  sets = paramDf[,list(set)][[1]]
  merge_one <- function(set,label,chr,stats,regions)
  {
    ourRegions = regions[[label]][[set]][[chr]]
    elementMetadata(ourRegions) = c(elementMetadata(ourRegions),stats)
    return(ourRegions)    
  }
  ourRegions = mclapply(1:nrow(paramDf),function(j){
    out = merge_one(sets[j],labels[j],chrs[j],merged_stats[[j]],regions)
    return(out)},mc.cores= 24)
  ourRegions = do.call(c,ourRegions)
  start(ourRegions) = pmin(start(ourRegions),ourRegions$newstart,na.rm=TRUE)
  end(ourRegions) = pmin(end(ourRegions),ourRegions$newend,na.rm=TRUE)
  p  = length(elementMetadata(ourRegions))
  elementMetadata(ourRegions) = elementMetadata(ourRegions)[1:(p-2)]
  return(ourRegions)
}

AY552regions = merge_all(merged_stats1,paramDf1[strand == "fwd"],regions)
AY553regions = merge_all(merged_stats2,paramDf2[strand == "fwd"],regions)
AY554regions = merge_all(merged_stats3,paramDf3[strand == "fwd"],regions)




save(list = "AY552regions",file = file.path(dr,"AY552_islandWsummary_lw3.RData"))
save(list = "AY553regions",file = file.path(dr,"AY553_islandWsummary_lw3.RData"))
save(list = "AY554regions",file = file.path(dr,"AY554_islandWsummary_lw3.RData"))



## summary_statistics <- function(reads)
## {
##   reads = data.table(as.data.frame(reads))
##   setkey(reads,strand)
##   fwd = na.omit(reads["+"])
##   bwd = na.omit(reads["-"])
##   f = nrow(fwd)
##   r = nrow(bwd)
  
##   f_pos = length(unique(fwd$start))
##   b_pos = length(unique(bwd$end))

##   max_fwd = NA
##   fwd_summit_pos = NA
  
##   if( f > 0){
##     fwd = GRanges(fwd$seqnames,ranges = IRanges(start = fwd$start,end=fwd$end),strand = "+")
##     fwd_cover = coverage(fwd)[[1]]
##     if(nrun(fwd_cover) > 1){
##       max_fwd = max(fwd_cover)
##       fwd_summit_pos = head(which(fwd_cover == max_fwd),n=1)
##     }
##   }

##   max_bwd = NA
##   bwd_summit_pos = NA
##   if(r > 0){  
##     bwd = GRanges(bwd$seqnames,ranges = IRanges(start = bwd$start,end=bwd$end),strand = "-")
##     bwd_cover = coverage(bwd)[[1]]
##     if(nrun(bwd_cover) > 1){
##       max_bwd = max(bwd_cover)
##       bwd_summit_pos = tail(which(bwd_cover == max_bwd),n=1)
##     }
##   }
##   depth = f + r
##   prob = f / depth
##   nrPos = f_pos + b_pos
  
##   diff = bwd_summit_pos - fwd_summit_pos    
##   out = c(f=f,r=r,prob=prob,depth=depth,f_uniq = f_pos,r_uniq = b_pos,n_uniq = nrPos,max_f = max_fwd,
##     max_r = max_bwd,max_f_pos = fwd_summit_pos,max_r_pos = bwd_summit_pos,diff = diff)
##   return(out)
    
## }

## both_idx = lapply(reg,function(x)as.character(which(x$label == "both")))

## sub_sep_reads_idx = mapply(function(x,y){
##   return(which(names(x) %in% y))},sep_reads,both_idx,SIMPLIFY=FALSE)

## sub_sep_reads = mapply(function(x,y){
##   return(x[y])},sep_reads,sub_sep_reads_idx,SIMPLIFY=FALSE)

## summary_stats = lapply(sub_sep_reads,function(x) data.table(t(mcmapply(summary_statistics,x,mc.cores=8))))

## save(file = file.path(datadir,"Ren2_summary_both.RData"),list = "summary_stats")
