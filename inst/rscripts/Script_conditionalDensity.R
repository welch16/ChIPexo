
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/chip.seq.set.from.pet.RData")
load("data/sample.summary.RData")

source("R/density_functions.R")
source("R/conditional_density_functions.R")

# param
figsdir = "inst/figs/cond_densities"
binsizes = c(300,500,750)
bins = lapply(binsizes,create.bins,genomeLength=seqlengths(exo[[1]]))
names(bins) = as.character(binsizes)


counts_strand <- function(bins,reads)
{
  counts_F = countOverlaps(bins,subset(reads,subset = strand(reads) == "+"))
  counts_R = countOverlaps(bins,subset(reads,subset = strand(reads) == "-"))
  return(DataFrame(counts_F,counts_R))
}


exo_strand_counts = mclapply(exo,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)
pet_strand_counts = mclapply(pet,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)
set_strand_counts = mclapply(set,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)

exo_strand_counts = unlist(exo_strand_counts)
pet_strand_counts = unlist(pet_strand_counts)
set_strand_counts = unlist(set_strand_counts)

exo_counts = mclapply(exo_strand_counts,function(x)do.call("+",x),mc.cores =6)
pet_counts = mclapply(pet_strand_counts,function(x)do.call("+",x),mc.cores =6)
set_counts = mclapply(set_strand_counts,function(x)do.call("+",x),mc.cores =6)

probs = c(0.9,0.95,0.975,.99)
exo_quantiles = mclapply(exo_counts,function(x)floor(quantile(x,probs)),mc.cores = 6)

pet_quantiles = mclapply(pet_counts,function(x)floor(quantile(x,probs)),mc.cores = 6)

set_quantiles = mclapply(set_counts,function(x)floor(quantile(x,probs)),mc.cores = 6)


fwd_strand_ratio <- function(f,r)return((f+1)/(f+r+2))

conditional_DF <- function(strand_counts,counts,quantile_vec)
{
  idx = lapply(quantile_vec,function(x)which(counts>x))
  ratio = lapply(idx,function(x)fwd_strand_ratio(strand_counts[x,1],strand_counts[x,2]))
  densities = lapply(ratio,density)
  out = lapply(densities,function(x)DataFrame(x$x,x$y))
  return(out)
}

exo_cond_densities = unlist(mcmapply(conditional_DF,exo_strand_counts,exo_counts,exo_quantiles,SIMPLIFY = FALSE,mc.cores = 6))

pet_cond_densities = unlist(mcmapply(conditional_DF,pet_strand_counts,pet_counts,pet_quantiles,SIMPLIFY = FALSE,mc.cores = 6))

set_cond_densities = unlist(mcmapply(conditional_DF,set_strand_counts,set_counts,set_quantiles,SIMPLIFY = FALSE,mc.cores = 6))

separate_name <- function(nm)
{
  split = strsplit(nm,"_")[[1]]
  edsn = split[1]
  split = strsplit(split[2],".",fixed=TRUE)[[1]]
  bin = split[2]
  prob = split[3]
  return(c(edsn,bin,prob))
}

assign_name <- function(dens_df,name)
{ 
  split = separate_name(name)
  dens_df$edsn = split[1]
  dens_df$bin = split[2]
  dens_df$prob = gsub("%","",split[3],fixed=TRUE)
  return(dens_df)
}

exo_cond_densities = mcmapply(assign_name,exo_cond_densities,names(exo_cond_densities),SIMPLIFY= FALSE,mc.cores =6)

pet_cond_densities = mcmapply(assign_name,pet_cond_densities,names(exo_cond_densities),SIMPLIFY= FALSE,mc.cores =6)

set_cond_densities = mcmapply(assign_name,set_cond_densities,names(exo_cond_densities),SIMPLIFY= FALSE,mc.cores =6)

save(list = c("exo_cond_densities","pet_cond_densities","set_cond_densities"),file = "data/cond_densities.RData")



exo_cond_densities = mclapply(exo_cond_densities,as.data.frame,mc.cores=6)
pet_cond_densities = mclapply(pet_cond_densities,as.data.frame,mc.cores=6)
set_cond_densities = mclapply(set_cond_densities,as.data.frame,mc.cores=6)


plot_cond_density <- function(df)
{
  p1 = ggplot(df,aes(x.x,x.y,colour = prob))+geom_line()+geom_vline(xintercept = .5,linetype="dashed")+
    facet_grid(edsn~bin,scales = "free")+
      theme(strip.text.y = element_text(angle = 0),legend.position ="bottom")+scale_x_continuous(limits = c(0,1))+xlab("")+ylab("")
  return(p1)
}

fix_set <- function(cond_densities)
{
  cond_density = do.call(rbind,cond_densities)
  cond_density$edsn = factor(cond_density$edsn)
  cond_density$bin  = factor(cond_density$bin )
  cond_density$prob = factor(cond_density$prob)
  return(cond_density)
}  

exo_cond_density = fix_set(exo_cond_densities)
pet_cond_density = fix_set(pet_cond_densities)
set_cond_density = fix_set(set_cond_densities)


plot_cond_density(set_cond_density)

pdf(file = file.path(figsdir,"chip_exo_new.pdf"),width = 5,height = 3*length(levels(exo_cond_density$edsn)))
p1 = plot_cond_density(exo_cond_density)
print(p1)
dev.off()

pdf(file = file.path(figsdir,"chip_pet_new.pdf"),width = 5,height = 3*length(levels(pet_cond_density$edsn)))
p2 = plot_cond_density(pet_cond_density)
print(p2)
dev.off()

pdf(file = file.path(figsdir,"chip_set_new.pdf"),width = 5,height = 3*length(levels(set_cond_density$edsn)))
p3 = plot_cond_density(set_cond_density)
print(p3)
dev.off()


# quantiles for table
# pet.quantiles = floor(do.call(rbind,lapply(pet,FUN = seq.quantile,binSize,probs)))
# exo.quantiles = floor(do.call(rbind,lapply(exo,FUN = seq.quantile,binSize,probs)))


## exo.depths = do.call(rbind,lapply(exo,FUN = length))
## pet.depths = do.call(rbind,lapply(pet,FUN = length))


# comparison exo vs pet of new data (edsn >1300)

## ip = c("Sig70","BetaPrimeFlag")
## rif = c("0 min","20 min")
## growth = "Aerobic"
## phase = "Exponential"
## j = 1
## st = list()

## L = length(ip)*length(rif)*length(growth)*length(phase)
## conditions =rep("",L)

## for(i in ip){
##   for(r in rif){
##     st[[j]] = resume.samples(ip = i,rif = r,growth = growth,phase = phase)
##     conditions[j] = paste0("ip",i,"_","rif",gsub(" ","",r))
##     j=j+1
##   }
## }

## for(i in 1:L){
##   tab = subset(sample.info,eval(parse(text = st[[i]])))
##   edsn = as.character(tab$edsn)
##   exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
##   exo.sets = lapply(exo.sets,function(y,exo)as(exo[[y]],"GRanges"),exo)
##   pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
##   pet.sets = lapply(pet.sets,function(y,pet)as(pet[[y]],"GRanges"),pet)
##   if(i == L){
##     pet.sets = list(pet.sets[[2]])
##   }
##   df_exo = quantile.multi.df(binsizes,probs,exo.sets,"exo",mc=8) 
##   df_pet = quantile.multi.df(binsizes,probs,pet.sets,"pet",mc=8)
##   df = rbind(df_exo,df_pet)
##   df$binSize = factor(df$binSize)
##   df$quantile = factor(df$quantile)
##   df$Rep = factor(df$Rep)
##   pdf(file = file.path(figsdir,paste0("Conditional_density_",conditions[i],".pdf")),width = 8,height =5)
##   for(prob in probs){
##     df1 = subset(df,quantile == prob )
##     p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = Rep))+geom_line()+facet_grid(type  ~ binSize,scales = "free")+
##       theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
##         scale_x_continuous(limits = c(0,1))+ggtitle(paste0("Quantile: ",prob*100, "%")) +
##         geom_vline(xintercept=.5,linetype = 2)
##     print(p)
##   }
##   dev.off() 
## }


# All chip exo forward strand ratio densities
exo_sample1 = names(exo)[grepl("13",names(exo))]
exo_sample2 = names(exo)[!grepl("13",names(exo))]


exo.sets1 = mclapply(exo_sample1,function(y,exo)as(exo[[y]],"GRanges"),exo,mc.cores =6)
exo.sets2 = mclapply(exo_sample2,function(y,exo)as(exo[[y]],"GRanges"),exo,mc.cores =6)

df_exo1 = quantile.multi.df(binsizes,probs,exo.sets1,"exo",mc=8) 
df_exo1$quantile = factor(df_exo1$quantile)
df_exo1$Rep =  factor(df_exo1$Rep,labels = exo_sample1)
names(df_exo1) = ifelse(names(df_exo1)=="Rep","sample",names(df_exo1))

df_exo2 = quantile.multi.df(binsizes,probs,exo.sets2,"exo",mc=8) 
df_exo2$quantile = factor(df_exo2$quantile)
df_exo2$Rep =  factor(df_exo2$Rep,labels = exo_sample2)
names(df_exo2) = ifelse(names(df_exo2)=="Rep","sample",names(df_exo2))
  
pdf(file = file.path(figsdir,"edsn1300_ChIP_exo_samples.pdf"),width = 8,height =20)
for(prob in probs){
  df1 = subset(df_exo1,quantile == prob )
  p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density))+geom_line()+facet_grid(sample  ~ binSize,scales = "free")+
    theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
      scale_x_continuous(limits = c(0,1))+ ggtitle(paste0("Quantile: ",prob*100, "%")) +
      geom_vline(xintercept=.5,linetype = 2)
  print(p)
}
dev.off()

pdf(file = file.path(figsdir,"edsn900_ChIP_exo_samples.pdf"),width = 8,height =20)
for(prob in probs){
  df1 = subset(df_exo2,quantile == prob )
  p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density))+geom_line()+facet_grid(sample  ~ binSize,scales = "free")+
    theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
      scale_x_continuous(limits = c(0,1))+ ggtitle(paste0("Quantile: ",prob*100, "%")) +
      geom_vline(xintercept=.5,linetype = 2)
  print(p)
}
dev.off()
  

