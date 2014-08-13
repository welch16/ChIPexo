
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/sample.summary.RData")

source("R/density_functions.R")
source("R/conditional_density_functions.R")

probs = c(.5,.75,.9,0.95,.975,.99,.999,0.9999,1)

binSize = 500

pet.quantiles = floor(do.call(rbind,lapply(pet,FUN = seq.quantile,binSize,probs)))
exo.quantiles = floor(do.call(rbind,lapply(exo,FUN = seq.quantile,binSize,probs)))

ip = c("Sig70","BetaPrimeFlag")
rif = c("0 min","20 min")
growth = "Aerobic"
phase = "Exponential"
j = 1
st = list()
for(i in ip){
  for(r in rif){
    st[[j]] = resume.samples(ip = i,rif = r,growth = growth,phase = phase)
    j=j+1
  }
}

exo.length = do.call(rbind,lapply(exo,FUN = length))
pet.length = do.call(rbind,lapply(pet,FUN = length))


quantile.df <- function(binSize,prob,type,Rep ,seqset)
{  
  bins = create.bins(binSize,seqlengths(seqset))
  bins$counts = countOverlaps(bins,seqset)
  seqquantile = quantile(bins$counts,prob)
  dens = density.reads.per.strand.ratio(subset(bins,subset = counts > seqquantile),
    seqset)
  df = data.frame("Fwd.Strand.Ratio"=dens$x,density =dens$y,binSize = binSize,quantile = prob,Rep=Rep,type = type)
  return(df)
}


quantile.multi.df <- function(binsizes,probs,exo.sets,pet.sets)
{
  
  all.binsizes = rep(binsizes,times = length(probs))
  all.prob = rep(probs,each =length(binsizes))

  exo.list = lapply(1:length(exo.sets),function(i,seqsets,all.binsizes,all.prob,type){
    df.list = mapply(quantile.df,all.binsizes,all.prob,MoreArgs = list(type,i,seqsets[[i]]),SIMPLIFY = FALSE)
    return(do.call(rbind,df.list))},exo.sets,all.binsizes,all.prob,"exo")

  pet.list = lapply(1:length(pet.sets),function(i,seqsets,all.binsizes,all.prob,type){
    df.list = mapply(quantile.df,all.binsizes,all.prob,MoreArgs = list(type,i,seqsets[[i]]),SIMPLIFY = FALSE)
    return(do.call(rbind,df.list))},pet.sets,all.binsizes,all.prob,"pet")

  df = rbind(do.call(rbind,exo.list),do.call(rbind,pet.list))
  
  df$binSize =factor(df$binSize)
  df$quantile = factor(df$quantile)
  df$Rep =factor(df$Rep)
  df$type = factor(df$type)

  return(df)
}

tab = subset(sample.info,eval(parse(text = st[[1]])))
edsn = as.character(tab$edsn)

exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
exo.sets = lapply(exo.sets,function(y,exo)as(exo[[y]],"GRanges"),exo)
pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
pet.sets = lapply(pet.sets,function(y,pet)as(pet[[y]],"GRanges"),pet)



binsizes = c(200,500,750)
probs = c(.3,0.5,.75,.9,.95)
df = quantile.multi.df(binsizes,probs,exo.sets,pet.sets)

pdf(width = 8,height =8)
df1 = subset(df,Rep==1)
p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = type))+geom_line()+facet_grid(quantile ~ binSize)
print(p)
df2 = subset(df,Rep==2)
p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = type))+geom_line()+facet_grid(quantile ~ binSize)
print(p)
dev.off()

