

rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(wavelets)
library(RColorBrewer)
library(scales)

load("data/chip.exo.RData")
load("data/sample.summary.RData")

source("R/depth_functions.R")

exo.edsn = c(1312,1315,1318,1321)
exo = exo[do.call(c,lapply(exo.edsn,function(x)grep(x,names(exo))))]

exo.edsn = exo.edsn[c(2,4)]


mc.cores = 8
exo = mclapply(exo,as.GRanges,mc.cores=mc.cores)


# this function calculates the number of reads in a distance-window for each unique fragment on set
d_ball_counts <- function(set,dist,mc.cores =8)
{
  posF = unique(start(subset(set,strand(set) == "+")))
  posR = unique(end(subset(set, strand(set) == "-")))
  pos = c(posF,posR)
  idx = 1:length(pos)  
  range = GRanges(seqnames= "U00096",range = IRanges(start = pos - dist,end = pos + dist),strand = "*")
  idx_split= split(idx,as.factor(sort(rep(1:80,length.out = length(pos))))) 
  partition = mclapply(idx_split,function(x,range)range[x],range,mc.cores =mc.cores)
  partition = mclapply(partition,function(x,set){
    counts = countOverlaps(x,set)
   x$counts = counts
   return(x)},set,mc.cores = mc.cores)
  return(unlist(GRangesList(partition))$counts)
}


ss = exo[[2]]
distances = c("25","50","100","250","500","1000")



d_counts = lapply(distances,function(x,set,width)d_ball_counts(set,as.numeric(x),ss,width)

## system.time(d_counts[["10"]] <- d_ball_counts(ss,10,width))
## system.time(d_counts[["25"]] <- d_ball_counts(ss,25,width))
## system.time(d_counts[["50"]] <- d_ball_counts(ss,50,width) )
## system.time(d_counts[["75"]] <- d_ball_counts(ss,75,width))
## system.time(d_counts[["100"]] <- d_ball_counts(ss,100,width))
## system.time(d_counts[["250"]] <- d_ball_counts(ss,250,width))
## system.time(d_counts[["500"]] <- d_ball_counts(ss,500,width))
## system.time(d_counts[["750"]] <- d_ball_counts(ss,750,width))
## system.time(d_counts[["1000"]] <- d_ball_counts(ss,1000,width))

## names(d_counts) = paste0("Ball_counts_",names(d_counts))
## save("d_counts",file = "Balls.RData")



pdf(file = "edsn1310_hexbinplot.pdf")
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(16)
dc = as.data.frame(d_counts)
for(i in 1:(length(distances)-1)){
  for(j in (i+1):length(distances)){
    d1 = paste0("X",distances[i])
    d2 = paste0("X",distances[j])
    p1 = ggplot(dc,aes(x=dc[[d1]],y=dc[[d2]]))+stat_binhex(bins = 100)+
      theme(legend.position = "bottom")+
      scale_x_log10( breaks=10^seq(0, 6, 2),labels=trans_format('log10', math_format(10^.x)) ) +
      scale_y_log10( breaks=10^seq(0, 6, 2),labels=trans_format('log10', math_format(10^.x)) ) +
      scale_fill_gradientn(limits = c(1,0.5*10^5),labels = trans_format('log10',math_format(10^.x)),colours=r, trans='log10')+
      xlab(d1)+ylab(d2)
    print(p1)
  }
}
dev.off() 

