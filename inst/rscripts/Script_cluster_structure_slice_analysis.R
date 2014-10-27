
rm(list = ls())
library(GenomicAlignments)
library(reshape2)
library(ggplot2)

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

alldata = mclapply(file.path(datadir,files),load_out,mc.cores =4)
baseDf = mcmapply(create_df,alldata,ip,rif,SIMPLIFY=FALSE,mc.cores=4)

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
  df = as.data.frame(do.call(rbind,df_list))
  df = melt(df)
  names(df)[6:7] = c("statistic","value")
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

levels(df$statistic)


figsdir = "inst/figs/summary_stats"
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
q2 = ggplot(subset(df,statistic == "diff" & label == "bwd" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Dark2")+ylab("Bwd. summit pos - Fwd summit pos")
print(q2)
q3 = ggplot(subset(df,statistic == "diff" & label == "fwd" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Dark2")+ylab("Bwd. summit pos - Fwd summit pos")
print(q3)
q4 = ggplot(subset(df,statistic == "diff" & seq == "exo" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Set1")+ylab("Bwd. summit pos - Fwd summit pos")
print(q4)
q4 = ggplot(subset(df,statistic == "diff" & seq == "pet" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Set1")+ylab("Bwd. summit pos - Fwd summit pos")
print(q4)
q6 = ggplot(subset(df,statistic == "diff" & seq == "set" ),aes(value,fill = interaction(seq,label)))+geom_density(alpha = I(.5))+facet_grid(rep~dataset)+theme(legend.position = "bottom")+scale_x_continuous(limits = c(-250,250))+geom_vline(xintercept = 0,linetype= "dashed")+scale_fill_brewer(palette = "Set1")+ylab("Bwd. summit pos - Fwd summit pos")
print(q6)
dev.off()

pdf(file = file.path(figsdir,"widthSlice.pdf"),width = 6,height =4)
z1 = ggplot(subset(df,statistic == "width" & seq == "exo" ),aes(value))+geom_freqpoly( aes(colour = label))+scale_x_continuous(limits = c(0,500))+theme(legend.position = "bottom")
print(z1)
dev.off()
