

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

exo.edsn = c(1312,1315,1318,1321)
pet.edsn = c(1397,1399,1401,1403)



exo = exo[do.call(c,lapply(exo.edsn,function(x)grep(x,names(exo))))]
pet = pet[do.call(c,lapply(pet.edsn,function(x)grep(x,names(pet))))]

# pick 20 min rif
exo.edsn = exo.edsn[c(2,4)]
pet.edsn = pet.edsn[c(2,4)]


## exo.edsn = exo.edsn[c(1,3)]
## pet.edsn = pet.edsn[c(1,3)]


exo = exo[do.call(c,lapply(exo.edsn,function(x)grep(x,names(exo))))]
pet = pet[do.call(c,lapply(pet.edsn,function(x)grep(x,names(pet))))]


mc.cores = 8
exo = mclapply(exo,as.GRanges,mc.cores=mc.cores)
pet = mclapply(pet,as.GRanges,mc.cores=mc.cores)
#set = mclapply(set,as.GRanges,mc.cores=mc.cores)


# using the second rep of chip exo sig70 - rif 20 because it has the highest
# depth
l = sapply(exo,length)
ss = exo[[which.max(l)]]
pos = start(ss)  

#U = sort(sample(1:seqlengths(ssF),2200))

distance = 500
idx = which(diff(pos) > distance)
length(idx)

## lowerBounds = U -500
## upperBounds = U + 500


lowerBounds = pmax(1,start(ss[idx])-distance)
upperBounds = pmin(seqlengths(ss[idx]),start(ss[idx]) + distance)


filterReads <- function(lb,ub,reads.pos)
  return(which(reads.pos >= lb & reads.pos <= ub))

exo.idx = list()
exo.idx[[1]] = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(exo[[1]])),SIMPLIFY = FALSE)
exo.idx[[2]] = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(exo[[2]])),SIMPLIFY = FALSE)

pet.idx = list()
pet.idx[[1]] = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(pet[[1]])),SIMPLIFY = FALSE)
pet.idx[[2]] = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(pet[[2]])),SIMPLIFY = FALSE)

## reads.idx.F = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(ssF)),SIMPLIFY = FALSE,
##   mc.cores = mc.cores)

## reads.idx.R1 = mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(end(ssR)),SIMPLIFY = FALSE,
##   mc.cores = mc.cores)

#snr <- function(reads_F,reads_R)return( (1 +length(reads_F))/(2 + length(reads_F) + length(reads_R)))
#SNR = do.call(c,mcmapply(FUN = snr,reads.idx.F,reads.idx.R,SIMPLIFY = FALSE,mc.cores=mc.cores))

getReads <- function(idx,reads,mc)mclapply(idx,function(x,reads)reads[x],reads,mc.cores =mc)

exo.reads = mapply(FUN = getReads,exo.idx,exo,MoreArgs = list(mc.cores),SIMPLIFY = FALSE)
pet.reads = mapply(FUN = getReads,pet.idx,pet,MoreArgs = list(mc.cores),SIMPLIFY = FALSE)

coverToVec <- function(lb,ub,cover)
{
  x = seq(lb,ub,by=1)
  if(nrun(cover)==1){
    y = rep(runValue(cover),length(x))
  }else{
    xf = cumsum(runLength(cover)[1:(nrun(cover)-1)])
    yf = runValue(cover)
    y = stepfun(xf,yf)(x)
  }
  return(y)
}


single_set_plot <- function(lb,ub,reads,main="")
{
  reads_F = subset(reads,strand(reads)=="+")
  reads_R = subset(reads,strand(reads)=="-")
  window_f= coverToVec(lb,ub,coverage(reads_F)[[1]])
  window_r= coverToVec(lb,ub,coverage(reads_R)[[1]])
  x = seq(lb,ub,by=1)
  xlim = c(lb,ub)
  ylim = c(-1,1) * max(max(window_f),max(window_r))
  plot(x =1 ,y=0,xlim = xlim,ylim = ylim,type = "n",main =main ,xlab = "",ylab = "" )  
  abline(h=0,col =  "black",lty =2)
  lines(c(lb,x,ub),c(0,window_f,0),col = "red")
  lines(c(lb,x,ub),c(0,-window_r,0),col = "blue")  
}

all_plots <- function(lb,ub,exo.reads1,exo.reads2,pet.reads1,pet.reads2)
{
  par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
  single_set_plot(lb,ub,exo.reads1,"exo1")
  single_set_plot(lb,ub,exo.reads2,"exo2")
  single_set_plot(lb,ub,pet.reads1,"pet1")
  single_set_plot(lb,ub,pet.reads2,"pet2")
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


pdf(file = "BetaPrimeFlag_0min.pdf")
for(i in 1:length(lowerBounds)){
  message(i)
  lb = lowerBounds[i]
  ub = upperBounds[i]
  all_plots(lb,ub,exo.reads[[1]][[i]],exo.reads[[2]][[i]],pet.reads[[1]][[i]],pet.reads[[2]][[i]])
}
dev.off()

pdf(file = "Rplots.pdf")
mcmapply(FUN = all_plots,lowerBounds[1:8],upperBounds[1:8],exo.reads[[1]][1:8],exo.reads[[2]][1:8],pet.reads[[1]][1:8],pet.reads[[2]][1:8],SIMPLIFY = FALSE);dev.off()         





## par(mfrow = c(2,2))
## x = seq(-5,5,length.out = 100)
## y = dunif(x)
## y = y / norm(as.matrix(y))
## plot(x,y,type = "l")
## yf = fft(y)
## plot(x, -log(Mod(yf)),type = "l")
## plot(x,Re(yf),type = "l")
## plot(x,Im(yf),type = "l")
## par(mfrow = c(2,2))
## x = seq(0,1,length.out = 100)
## y = dunif(x)
## y = y / norm(as.matrix(y))
## plot(x,y,type = "l")
## yf = fft(y)
## plot(x, -log(Mod(yf)),type = "l")
## plot(x,Re(yf),type = "l")
## plot(x,Im(yf),type = "l")
## dev.off()




