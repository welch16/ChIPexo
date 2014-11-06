

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
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
figsdir = "inst/figs/peak_plots"
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
set.sets = lapply(edsn,function(x)set[do.call(c,lapply(x,function(y)grep(y,names(set))))])


exo.sets = lapply(exo.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)
pet.sets = lapply(pet.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)
set.sets = lapply(set.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)

exo.depths = sapply(exo.sets,function(x)c(rep1=length(x[[1]]),rep2=length(x[[2]])))
pet.depths = sapply(pet.sets,function(x)c(rep1=length(x[[1]]),rep2=length(x[[2]])))
set.depths = sapply(set.sets,function(x)c(rep1=length(x[[1]]),rep2=length(x[[2]])))

max.depth.idx = ifelse(exo.depths[1,] > exo.depths[2,],1,2)
exo.sets.use = mapply(function(x,y) x[[y]],exo.sets,max.depth.idx,SIMPLIFY =FALSE)

fwd_cover = mclapply(exo.sets.use,function(x)coverage(subset(x,strand(x)=="+")))
bwd_cover = mclapply(exo.sets.use,function(x)coverage(subset(x,strand(x)=="-")))



## i=4
## v = summary(runValue(fwd_cover[[i]])[[1]])
## w = summary(runValue(bwd_cover[[i]])[[1]])
## sapply(v,function(x)length(slice(fwd_cover[[i]],upper = x)[[1]]))
## sapply(w,function(x)length(slice(bwd_cover[[i]],upper = x)[[1]]))

# note: the min numer of regions is given by the mean..
fwd_regions = mclapply(fwd_cover,function(x)slice(x,lower = 5,rangesOnly=TRUE)[[1]],mc.cores=4)  
bwd_regions = mclapply(bwd_cover,function(x)slice(x,lower = 5,rangesOnly=TRUE)[[1]],mc.cores=4)


# convert to GRanges
fwd_regions = mclapply(fwd_regions,function(x)GRanges(seqnames ="U00096",ranges =x ,strand = "*"),mc.cores=4)
bwd_regions = mclapply(bwd_regions,function(x)GRanges(seqnames ="U00096",ranges =x ,strand = "*"),mc.cores=4)

minW = 0
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

# We are gonna try first considering all regions

ourRegions = mapply(create_regions,fwd_regions,bwd_regions,MoreArgs = list(minW=0))

# First plot for analysis, compare the behaviour of the width of the regions
ip = rep(c("Sig70","BetaPrimeF"),2)
rif = rep(paste(c(0,20),"min"),each=2)

df = do.call(rbind,mapply( function(x,y,z){
  df = data.frame(width = width(x))
  df$label = x$label
  df$ip = y
  df$rif = z
  return(df)},ourRegions,ip,rif,SIMPLIFY=FALSE)
)

df$treatment = paste0(df$ip,"-",df$rif)

# I think the shape of the plot is happening the way I am building the
# ranges, since I am using the union operator for one case and not the
# other then

pdf(file = file.path(figsdir,"widthPrev.pdf"))
p1 = ggplot(df,aes(treatment,width,colour = label))+geom_boxplot() +facet_grid(.~treatment,drop=TRUE,scales="free_x")+ scale_y_log10()
print(p1)
dev.off()

separated_reads<- function(reads,regions)
{
  overlaps = findOverlaps(regions,reads)
  sep_index = split(overlaps,factor(queryHits(overlaps)))
  separated_reads <- mclapply(sep_index,
    function(x,reads)reads[subjectHits(x)],reads,mc.cores = 8)
  return(separated_reads)
}

datadir = "data"
rif = c(0,0,20,20)
files = paste0(ip,"_",rif,"min")
ext= 100

save(list = "ourRegions", file = file.path(datadir,"regions_slice.RData"))


for(i in 1:4){
  print(i)
  regions = ourRegions[[i]]
  print("Separating exo reads into regions")
  exo_sep_reads1 = separated_reads(exo.sets[[i]][[1]],regions)
  exo_sep_reads2 = separated_reads(exo.sets[[i]][[2]],regions)
  print("Separating pet reads into regions")
  pet_sep_reads1 = separated_reads(pet.sets[[i]][[1]],regions)
  pet_sep_reads2 = separated_reads(pet.sets[[i]][[2]],regions)
  print("Separating set reads into regions")
  set_sep_reads1 = separated_reads(set.sets[[i]][[1]],regions)
  set_sep_reads2 = separated_reads(set.sets[[i]][[2]],regions)

  # note the list have different lengths
  
  save(list = c("exo_sep_reads1","exo_sep_reads2","pet_sep_reads1",
         "pet_sep_reads2","set_sep_reads1","set_sep_reads2"),
       file = file.path(datadir,paste0(files[i],"_sep_reads.RData")))
  
  lowerBounds = start(regions)
  upperBounds = end(regions)

  names(upperBounds) = as.character(1:length(regions))
  names(lowerBounds) = as.character(1:length(regions))
  
  pdf(file = file.path(figsdir,paste0(files[i],"_sliceRegions.pdf")),width =12,height =6)
  lapply(as.character(1:length(regions)),function(j){
    all_plots2(lowerBounds[j],upperBounds[j],
       exo_sep_reads1[[j]],exo_sep_reads2[[j]],
       pet_sep_reads1[[j]],pet_sep_reads2[[j]],
       set_sep_reads1[[j]],set_sep_reads2[[j]],ext)
  })
  dev.off()

}

