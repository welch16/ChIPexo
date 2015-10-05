
rm(list = ls())
library(GenomicAlignments)
library(reshape2)
library(ggplot2)
library(data.table)

datadir = "data"
files = c("outputSig70_0min_sep_reads.RData",
  "outputSig70_20min_sep_reads.RData",
  "outputBetaPrimeF_0min_sep_reads.RData",
  "outputBetaPrimeF_20min_sep_reads.RData")
ip = rep(c("Sig70","BetaPrimeF"),each=2)
rif = rep(c("0min","20min"),2)

load_out <- function(file)
{
  load(file)
  l1 = list(regions,exo1,exo2,pet1,pet2,set1,set2)
  names(l1) = c("regions","exo1","exo2","pet1","pet2","set1","set2")
  return(l1)
}



create_df <- function(alldata,ip,rif)
{
  regions = alldata[[1]]
  n = length(regions)
  df = elementMetadata(regions)
  df$width = width(regions)
  df$ip = ip
  df$rif = rif  
  return(df)
}

readfiles = c("Sig70_0min_sep_reads.RData",
  "Sig70_20min_sep_reads.RData",
  "BetaPrimeF_0min_sep_reads.RData",
  "BetaPrimeF_20min_sep_reads.RData")


alldata = mclapply(file.path(datadir,files),load_out,mc.cores =4)
baseDf = mcmapply(create_df,alldata,ip,rif,SIMPLIFY=FALSE,mc.cores=4)
load(file = file.path(datadir,"position_numer.RData"))
names(alldata) = readfiles
names(baseDf) =readfiles
names(posList) =readfiles

add_positions <- function(exp_data,positions)
{
  exp_data$exo1[["nrPos"]] = positions[[1]]
  exp_data$exo2[["nrPos"]] = positions[[2]]
  exp_data$exo1$readsPos_ratio = ifelse(positions[[1]]==0,NA,(exp_data$exo1[["f"]] + exp_data$exo1[["r"]]) / positions[[1]])
  exp_data$exo2$readsPos_ratio = ifelse(positions[[2]]==0,NA,(exp_data$exo2[["f"]] + exp_data$exo2[["r"]]) / positions[[2]])
  return(exp_data)
}


alldata = mcmapply(add_positions,alldata,posList,SIMPLIFY=FALSE,mc.cores =4)

merge_df <- function(baseDf,alldata)
{  
  df = DataFrame(chrID = seqnames(alldata$regions))
  df$start = start(alldata$regions)
  df$end = end(alldata$regions)
  df = cbind(df,baseDf)
  alldata = alldata[-1]
  alldata = mapply(function(x,y){
    x$seq = y
    return(x)},alldata,names(alldata))
  df_list = lapply(alldata,function(x,df)cbind(df,x),df)
  df_list = lapply(df_list,function(x)melt(as.data.frame(x)))
  df = do.call(rbind,df_list)
  names(df)[6:7] = c("statistic","value")
  df$statistic = as.character(df$statistic)
  return(df)
}


df_list = mcmapply(merge_df,baseDf,alldata,SIMPLIFY = FALSE,mc.cores = 4)

df = do.call(rbind,df_list)
df$rep = factor(as.numeric(substring(df$seq,4,4)))
df$seq = factor(substring(df$seq,1,3))
df$label = factor(df$label)
df$dataset = factor(paste0(df$ip,":",df$rif))
df$ip = factor(df$ip)
df$rif = factor(df$rif)


figsdir = "inst/figs/summaryStats"
pdf(file = file.path(figsdir,"FwdStrandRatioSlice.pdf"),height = 6,width =12)
p0 = ggplot(subset(df,statistic == "prob"),aes(seq,value))+geom_boxplot(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+ylab("Fwd. Strand ratio")
print(p0)
p1 = ggplot(subset(df,statistic == "prob"),aes(seq,value,fill = label))+geom_boxplot(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+ylab("Fwd. Strand ratio")
print(p1)
p2 = ggplot(subset(df,statistic == "prob"),aes(label,value,fill = seq))+geom_boxplot(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Dark2")+ylab("Fwd. Strand ratio")
print(p2)
dev.off()

pdf(file = file.path(figsdir,"SummitPosDiffSlice.pdf"),height = 6,width =12)
q1 = ggplot(subset(df,statistic == "diff" & label == "both" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Dark2")+ylab("Bwd. summit pos - Fwd summit pos")
print(q1)
q4 = ggplot(subset(df,statistic == "diff" & seq == "exo" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Set1")+ylab("Bwd. summit pos - Fwd summit pos")
print(q4)
dev.off()

pdf(file = file.path(figsdir,"widthSlice.pdf"),width = 6,height =4)
z1 = ggplot(subset(df,statistic == "width" & seq == "exo" ),aes(value))+geom_freqpoly( aes(colour = label))+scale_x_continuous(limits = c(0,500))+theme(legend.position = "bottom")
print(z1)
dev.off()

d = levels(df$dataset)

df = data.table(df)
setkey(df,label,ip,rif,seq,statistic,rep,dataset)

build_hexbin <- function(dataset,df)
{
  df1 = df[dataset == dataset & seq == "exo"]
  width_df = df1[statistic == "width"]
  ratio = df1[statistic == "readsPos_ratio"]$value
  prob = df1[statistic == "prob"]$value
  depth = df1[statistic == "f"]$value + df1[statistic =="r"]$value
  summitDiff = df1[statistic == "diff"]$value
  cols = c("label","dataset","rep","width")
  df1 = width_df[,list(label,dataset,rep,value)]
  setnames(df1,colnames(df1),cols)
  df1$posreads_ratio = ratio
  df1$depth = depth
  group = factor(colSums(do.call(rbind,lapply(seq(0.2,1,by=.2),function(x,prob)prob <= x,prob))))
  group = mapvalues(group,from=as.character(1:5),to=paste0(seq(0.8,0,by=-.2),"-",seq(1,.2,by=-.2)))
  df1$prob = prob
  df1$group = group
  df1$summitDiff = summitDiff   
  return(df1)
}

hb_df = mclapply(d,build_hexbin,df,mc.cores=4)

names(hb_df) = d
#hb_df = do.call(rbind,hb_df)

library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(RColorBrewer)
library(scales)

rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(16)

hb = do.call(rbind,hb_df)

maxwidth = 1000
maxratio = 100
nrbins = 60
readsLen = 51


q1 = ggplot(subset(hb,rep == 1 & label == "both") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset,nrow=2)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 1 (both)")+ylab("Nr. reads / Nr. positions")
print(q1)
q2 = ggplot(subset(hb,rep == 1 & label == "fwd") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 1 (fwd)")+geom_vline(xintercept = readsLen,linetype = "dashed")+ylab("Nr. reads / Nr. positions")
print(q2)
q3 = ggplot(subset(hb,rep == 1 & label == "bwd") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 1 (bwd)")+geom_vline(xintercept = readsLen,linetype = "dashed")+ylab("Nr. reads / Nr. positions")
print(q3)
q1 = ggplot(subset(hb,rep == 2 & label == "both") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset,nrow=2)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 2 (both)")+ylab("Nr. reads / Nr. positions")
print(q1)
q2 = ggplot(subset(hb,rep == 2 & label == "fwd") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 2 (fwd)")+geom_vline(xintercept = readsLen,linetype = "dashed")+ylab("Nr. reads / Nr. positions")
print(q2)
q3 = ggplot(subset(hb,rep == 2 & label == "bwd") ,aes(width,posreads_ratio))+stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+
  facet_wrap(~dataset)+scale_x_continuous(limits= c(0,maxwidth))+scale_y_continuous(limits = c(0,maxratio))+
  ggtitle("exo - rep 2 (bwd)")+geom_vline(xintercept = readsLen,linetype = "dashed")+ylab("Nr. reads / Nr. positions")
print(q3)



pdf(file = file.path(figsdir,"summaryhexbin.pdf"))
nrbins = 60
q5 = ggplot(hb[label=="both"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+ scale_fill_gradientn(colours=r,
  trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+    
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits  = 250* c(-1,1))+ggtitle("both")
print(q5)
q6 = ggplot(hb[label=="both"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+ scale_fill_gradientn(colours=r,
  trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+    
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits  = 250* c(-1,1))+ggtitle("both")
print(q6)

q5 =  ggplot(hb[label=="fwd"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+ theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size =   .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits =
  250* c(-1,1))+ggtitle("fwd")
print(q5)
q6 = ggplot(hb[label=="fwd"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+ theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,150))+ scale_y_continuous(limits =
  250* c(-1,1))+ggtitle("fwd")
print(q6)

q5 =  ggplot(hb[label=="bwd"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+ theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size =   .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits =
  250* c(-1,1))+ggtitle("bwd")
print(q5)
q6 = ggplot(hb[label=="bwd"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+ theme(legend.position = "bottom")+
  scale_fill_gradientn(colours=r, trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,150))+ scale_y_continuous(limits =
  250* c(-1,1))+ggtitle("bwd")
print(q6)


dev.off()



q7 = ggplot(hb[label!="both"],aes(depth,summitDiff))+stat_binhex(bins=nrbins)+
  stat_binhex(bins = nrbins)+
  theme(legend.position = "bottom")+ scale_fill_gradientn(colours=r,
  trans='log10')+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+    
  geom_abline(intercept =51,slope = 0,linetype = 3,size =.8)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =.8)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits  = 250* c(-1,1))
print(q7)



q7 = ggplot(hb[label!="both"],aes(depth,summitDiff,colour = group))+geom_point(alpha = I(.15),shape = 0)+
  theme(legend.position = "bottom")+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+       
  geom_abline(intercept =51,slope = 0,linetype = 2,size =1)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =1)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits  = 250* c(-1,1))+scale_color_brewer(palette="Set1")
q7

q8 = ggplot(hb[label =="fwd"],aes(depth,summitDiff,colour = group))+geom_point(alpha = I(.15),shape = 0)+
  theme(legend.position = "bottom")+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+       
  geom_abline(intercept =51,slope = 0,linetype = 2,size =1)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =1)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits  = 250* c(-1,1))+scale_color_brewer(palette="Set1")
q8

q9 = ggplot(hb[label=="bwd"],aes(depth,summitDiff,colour = group))+geom_point(alpha = I(.15),shape = 0)+
  theme(legend.position = "bottom")+facet_wrap(~dataset)+
  geom_smooth(se = FALSE,size = 1,aes(group=group,colour = group),method="loess")+scale_color_brewer(palette="Set1")+
  geom_smooth(se = FALSE,size = .8,aes(group=1),colour=I("black"),linetype=2,method="loess")+       
  geom_abline(intercept =51,slope = 0,linetype = 2,size =1)+
  geom_abline(intercept =0,slope = 0,linetype = 3,size =1)+
  scale_x_continuous(limits = c(0,500))+
  scale_y_continuous(limits  = 250* c(-1,1))+scale_color_brewer(palette="Set1")
q9




## stats = unique(df$statistic)
## df1 = subset(df,statistic == "f")
## df1$statistic = "depth"
## df1$value = df1$value + subset(df,statistic=="r")$value

## df = rbind(df,df1)
## stats = unique(df$statistic)


## df1 = subset(df,seq =="exo")

## depth1 = subset(df1,statistic == "depth" & rep == 1)



## levels(df1$dataset)





## load(file.path(datadir,readfiles[1])) # exo, pet and set: rep1 and rep2

## peak =2317
## exo1 = exo_sep_reads1[[peak]]
## exo2 = exo_sep_reads2[[peak]]
## pet1 = pet_sep_reads1[[peak]]
## pet2 = pet_sep_reads2[[peak]]
## set1 = set_sep_reads1[[peak]]
## set2 = set_sep_reads2[[peak]]

## ext= 100
## alldata[[1]]$regions[peak]


## base = baseDf[[1]][peak,]

## stats = as.data.frame(do.call(rbind,lapply(alldata[[1]][-1],function(x)x[peak,])))
## rownames(stats) = c("exo1","exo2","pet1","pet2","set1","set2")


