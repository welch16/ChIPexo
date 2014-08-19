
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

# param
binsizes = c(200,500,750,1000)
probs = c(0,0.5,.75,0.9,0.95,0.975)
figsdir = "inst/figs/cond_densities"


# quantiles for table
# pet.quantiles = floor(do.call(rbind,lapply(pet,FUN = seq.quantile,binSize,probs)))
# exo.quantiles = floor(do.call(rbind,lapply(exo,FUN = seq.quantile,binSize,probs)))


## exo.depths = do.call(rbind,lapply(exo,FUN = length))
## pet.depths = do.call(rbind,lapply(pet,FUN = length))


# comparison exo vs pet of new data (edsn >1300)
ip = c("Sig70","BetaPrimeFlag")
rif = c("0 min","20 min")
growth = "Aerobic"
phase = "Exponential"
j = 1
st = list()

L = length(ip)*length(rif)*length(growth)*length(phase)
conditions =rep("",L)

for(i in ip){
  for(r in rif){
    st[[j]] = resume.samples(ip = i,rif = r,growth = growth,phase = phase)
    conditions[j] = paste0("ip",i,"_","rif",gsub(" ","",r))
    j=j+1
  }
}

for(i in 1:L){
  tab = subset(sample.info,eval(parse(text = st[[i]])))
  edsn = as.character(tab$edsn)
  exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
  exo.sets = lapply(exo.sets,function(y,exo)as(exo[[y]],"GRanges"),exo)
  pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
  pet.sets = lapply(pet.sets,function(y,pet)as(pet[[y]],"GRanges"),pet)
  if(i == L){
    pet.sets = list(pet.sets[[2]])
  }
  df_exo = quantile.multi.df(binsizes,probs,exo.sets,"exo",mc=8) 
  df_pet = quantile.multi.df(binsizes,probs,pet.sets,"pet",mc=8)
  df = rbind(df_exo,df_pet)
  df$binSize = factor(df$binSize)
  df$quantile = factor(df$quantile)
  df$Rep = factor(df$Rep)
  pdf(file = file.path(figsdir,paste0("Conditional_density_",conditions[i],".pdf")),width = 8,height =5)
  for(prob in probs){
    df1 = subset(df,quantile == prob )
    p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = Rep))+geom_line()+facet_grid(type  ~ binSize,scales = "free")+
      theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
        scale_x_continuous(limits = c(0,1))+ggtitle(paste0("Quantile: ",prob*100, "%")) +
        geom_vline(xintercept=.5,linetype = 2)
    print(p)
  }
  dev.off() 
}


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
  p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = quantile))+geom_line()+facet_grid(sample  ~ binSize,scales = "free")+
    theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
      scale_x_continuous(limits = c(0,1))+ ggtitle(paste0("Quantile: ",prob*100, "%")) +
      geom_vline(xintercept=.5,linetype = 2)
  print(p)
}
dev.off()

  
## pdf(file = "newdata_ChIP_exo_samples2.pdf",width = 8,height =20)
## df1 = df_exo_exp
## p <- ggplot(df1,aes(x=Fwd.Strand.Ratio,y=density,colour = quantile))+geom_line()+facet_grid(sample  ~ binSize,scales = "free") +theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
##     scale_x_continuous(limits = c(0,1))+ geom_vline(xintercept=.5,linetype = 2)
## print(p)
## dev.off()






