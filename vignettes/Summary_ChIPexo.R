## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----loaddepth,include=FALSE,echo=FALSE,eval=TRUE---------------------------------------
  library(knitr)
  library(reshape2)
  load("../data/depth.RData")
  exo.depth = melt(exo.depth)
  exo.depth = exo.depth[,c(2,1)]
  names(exo.depth) = c("Data set","Depth")
  exo.depth[,1]= sapply(exo.depth[,1],function(x)gsub("_042814","",x))
  pet.depth = melt(pet.depth)
  pet.depth = pet.depth[,c(2,1)]
  pet.depth[,1]= sapply(pet.depth[,1],function(x)gsub("_042814","",x))
  names(pet.depth) = c("Data set","Depth")

## ----include=TRUE,echo=FALSE,eval=TRUE,results="asis"-----------------------------------
 kable(exo.depth,format = 'latex')

## ----include=TRUE,echo=FALSE,eval=TRUE,results="asis"-----------------------------------
  kable(pet.depth,format ='latex')

## ----load_cross_data,include=FALSE,echo=FALSE,eval=TRUE---------------------------------
  library(ggplot2)
  datadir = "../data"
  load(file.path(datadir,"cross_corr.RData"))

## ----chipseqSET,include=FALSE,echo=FALSE,eval=TRUE--------------------------------------
  df = set_df_list[[3]]       
  ggplot(df,aes(shift,crossCorr,colour = method))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1/3,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2)),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")
  ggsave("setCC1.pdf",width = 12,height = 5)

## ----chipseqPET,include=FALSE,echo=FALSE,eval=TRUE--------------------------------------
  df = pet_df_list[[5]]       
  ggplot(df,aes(shift,crossCorr,colour = method))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1/3,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2)),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")
  ggsave("setCC2.pdf",width = 12,height = 5)

## ----chipseqEXO,include=FALSE,echo=FALSE,eval=TRUE--------------------------------------
  df = exo_df_list[[8]]       
  ggplot(df,aes(shift,crossCorr,colour = method))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1/3,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2)),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")
  ggsave("setCC3.pdf",width = 12,height = 5)

## ----loadensities,include=FALSE,echo=FALSE,eval=TRUE------------------------------------
load(file.path(datadir,"densities.RData")) # exo_density, pet_density,set_density
load(file.path(datadir,"sample.summary.RData")) # sample.info
source("../R/density_functions.R")
row= subset(sample.info,eval(parse(text=resume.samples(edsn=1319))))
rows = subset(sample.info,(eval(parse(text = resume.samples(ip=row$ip,rif=row$rif)))))
rm(row)
idd = do.call(c,lapply(rows$edsn,function(x)grep(x,names(exo_density))))
rm(rows)

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
  df = do.call(rbind, lapply(exo_density[idd],as.data.frame))
  idd
  ggplot(df,aes(x.x,x.y,colour = bin))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2),angle=90),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")+facet_grid(.~edsn)+scale_x_continuous(limits = c(0,1))+xlab("fwd. strand ratio")+ylab("density")
  ggsave("exoDens.pdf",width = 6,height = 6)
  rm(idd)

## ----chipseq_dens,include=FALSE,echo=FALSE,eval=TRUE------------------------------------
row= subset(sample.info,eval(parse(text=resume.samples(edsn=1400))))
idd = do.call(c,lapply(row$edsn,function(x)grep(x,names(pet_density))))

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
  df1 = do.call(rbind,lapply(pet_density[idd],as.data.frame))
  df1$seq = "PET"
  df2 = do.call(rbind,lapply(set_density[idd],as.data.frame))
  df2$seq = "SET"
  df = rbind(df1,df2)
  df$seq = factor(df$seq)
  ggplot(df,aes(x.x,x.y,colour = bin))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2),angle=90),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")+facet_grid(edsn~seq)+scale_x_continuous(limits = c(0,1))+xlab("fwd. strand ratio")+ylab("density")
  ggsave("seqDens.pdf",width = 6,height = 6)  

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
load(file = "../data/cond_densities.RData")
fix_set <- function(cond_densities)
{
  cond_density = do.call(rbind,cond_densities)
  cond_density$edsn = factor(cond_density$edsn)
  cond_density$bin  = factor(cond_density$bin )
  cond_density$prob = factor(cond_density$prob)
  return(cond_density)
}
exo_cond_densities = lapply(exo_cond_densities,as.data.frame)
exo_cond_density = fix_set(exo_cond_densities)

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
  df = subset(exo_cond_density,edsn == "edsn1310")
  ggplot(df,aes(x.x,x.y,colour = prob))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2),angle=90),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")+facet_grid(edsn ~ bin)+scale_x_continuous(limits = c(0,1))+xlab("fwd. strand ratio")+ylab("density") 
  ggsave("exo1310.pdf",width = 12,height = 6)  

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
  df = subset(exo_cond_density,edsn == "edsn1319")
  ggplot(df,aes(x.x,x.y,colour = prob))+geom_line()+theme(legend.position = "bottom",aspect.ratio = 1,axis.title = element_text(size = rel(2.5)),axis.text = element_text(size = rel(2),angle=90),legend.text = element_text(size = rel(2)),legend.title = element_text(size = rel(2)))+ scale_color_brewer(palette ="Dark2")+facet_grid(edsn ~ bin)+scale_x_continuous(limits = c(0,1))+xlab("fwd. strand ratio")+ylab("density") 
  ggsave("exo1319.pdf",width = 12,height = 6)  

