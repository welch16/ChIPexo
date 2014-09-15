

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

exo.sets = lapply(exo.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)
pet.sets = lapply(pet.sets,function(x,mc)mclapply(x,as.GRanges,mc.cores=mc),mc.cores)


indexes = mapply(FUN = filterSets,exo.sets,pet.sets,MoreArgs = list(distance),SIMPLIFY = FALSE)
names(indexes) = conditions

save(list = "indexes",file = "data/indexes_regions_1.RData")

filenames = c("Sig70_rif0_peak_plots.pdf","BetaPrimeFlag_rif0_peak_plots.pdf",
  "Sig70_rif20_peak_plots.pdf","BetaPrimeFlag_rif20_peak_plots.pdf")

# plots
for(i in 1:4)
{  
  exo.reads = mapply(FUN = getReads,indexes[[i]]$exo,exo.sets[[i]],MoreArgs = list(mc.cores),SIMPLIFY = FALSE)
  pet.reads = mapply(FUN = getReads,indexes[[i]]$pet,pet.sets[[i]],MoreArgs = list(mc.cores),SIMPLIFY = FALSE)

  j = which.max(sapply(exo.sets[[i]],length))
  ss = exo.sets[[i]][[j]]

  pos = start(ss)
  idx = which(diff(pos) > distance)
  
  lowerBounds = pmax(1,start(ss[idx])-distance)
  upperBounds = pmin(seqlengths(ss[idx]),start(ss[idx]) + distance)
  
  pdf(file = file.path(figsdir,filenames[i]))
  for(k in 1:length(lowerBounds)){
    # message(k)
    lb = lowerBounds[k]
    ub = upperBounds[k]
    all_plots(lb,ub,exo.reads[[1]][[k]],exo.reads[[2]][[k]],pet.reads[[1]][[k]],pet.reads[[2]][[k]])
  }
  dev.off()

  write.csv(data.frame(lower = lowerBounds,upper = upperBounds),file = gsub(".pdf",".csv",filenames[i]),row.names = FALSE)
}



## i=113
## lb = lowerBounds[i]p
## ub = upperBounds[i]
## cover_f = cover_F[[i]]
## cover_r = cover_R[[i]]
## window_f = as.matrix(coverToVec(lb,ub,cover_f))
## dwt_f =  dwt(window_f,"haar",boundary ="reflection")
## window_r = as.matrix(coverToVec(lb,ub,cover_r))
## dwt_r =  dwt(window_r,"haar",boundary ="reflection")
## pdf(height = 3.5)
## plot.dwt.multiple(list(dwt_f,dwt_r),levels =3)
## for(i in 1:length(dwt_f@V)){plot(dwt_f@V[[i]] - dwt_r@V[[i]],type = "l");abline(h=0,col = "red",lty=2)
##                             plot(dwt_f@W[[i]] - dwt_r@W[[i]],type = "l");abline(h=0,col = "red",lty=2)}
## dev.off()


## pdf(file = "Rplots.pdf")
## mcmapply(FUN = all_plots,lowerBounds[1:8],upperBounds[1:8],exo.reads[[1]][1:8],exo.reads[[2]][1:8],pet.reads[[1]][1:8],pet.reads[[2]][1:8],SIMPLIFY = FALSE);dev.off()         





