
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(RColorBrewer)
library(scales)

load("data/sample.summary.RData")
load("data/summary_stats_labeled_regions_withSet.RData")

source("R/depth_functions.R")
source("R/density_functions.R")
source("R/structure_functions.R")

# Parameters
mc.cores = 8
distance = 500
figsdir = "inst/figs/summary_stats"


# Set conditions
conditions = rep("",4)
conditions[1] = resume.samples(ip="Sig70",rif="0 min")
conditions[2] = resume.samples(ip="BetaPrimeFlag",rif="0 min")
conditions[3] = resume.samples(ip="Sig70",rif="20 min")
conditions[4] = resume.samples(ip="BetaPrimeFlag",rif="20 min")

ip = rep(c("Sig70","BetaPF"),2)
rif = rep(c(0,20),each = 2)

metadata_to_df <- function(summary_table,ip,rif)
{
  table = elementMetadata(summary_table) 
  df = data.frame(ip = ip,rif = rif,start = rep(start(summary_table),4),end = rep(end(summary_table),4),
    isPeak = ifelse(rep(table[["isPeak"]],4)==1,"yes","no")   )
  df$seq = rep(c("exo","set"),each = length(summary_table)*2)
  df$rep = rep(1:2,each = length(summary_table))  
  df$countsScaled = c(table[["exo1scaled"]],table[["exo2scaled"]],table[["set1scaled"]],table[["set2scaled"]])
  df$fwdStrandRatio = c(table[["exo1_fwdStrandRatio"]],table[["exo2_fwdStrandRatio"]],
    table[["set1_fwdStrandRatio"]],table[["set2_fwdStrandRatio"]])
  df$summitPosDiff = c(table[["exo1_summit_pos_diff"]],table[["exo2_summit_pos_diff"]],
    table[["set1_summit_pos_diff"]],table[["set2_summit_pos_diff"]])
  return(df)
}

df_list= mcmapply(FUN = metadata_to_df,summary_stats_gr,ip,rif,SIMPLIFY = FALSE,mc.cores =4)
df = do.call(rbind,df_list)
df$rif = factor(df$rif,levels = c(0,20))
df$isPeak = factor(df$isPeak,levels = c("yes","no"))
df$experiment = paste0(df$ip,"-",df$rif,"min")
df$experiment = factor(df$experiment)
df$seq = factor(df$seq)
df$rep = factor(df$rep)


pdf(file = file.path(figsdir,"Fwd_strand_Ratio_boxplot.pdf"))
p1 = ggplot(df,aes(isPeak,fwdStrandRatio,colour = experiment))+geom_boxplot()+facet_grid(rep~seq)+ theme(legend.position = "bottom")+ylab("Forward strand ratio")
print(p1)
dev.off()


pdf(file = file.path(figsdir,"Strand_summit_position_boxplot.pdf"))
p2 = ggplot(df,aes(isPeak,summitPosDiff,colour = experiment))+geom_boxplot()+facet_grid(rep~seq)+ theme(legend.position = "bottom")+ylab("Difference in strand summit position")
print(p2)
dev.off()

pdf(file = file.path(figsdir,"sqrtDepthScaled_boxplot.pdf"))
p3 = ggplot(df,aes(isPeak,sqrt(countsScaled),colour = experiment))+geom_boxplot()+facet_grid(rep~seq)+theme(legend.position = "bottom")+ylab("Squared root of scaled counts")
print(p3)
dev.off()


df1 = subset(df,ip == 'BetaPF' &rif == 0 & rep == 1 & seq == "exo" &  isPeak == "yes")
which(df1$summitPosDiff > 490)
df1[23,]

which(df1$summitPosDiff < -400)
df1[c(75,83,139),]


df2 = subset(df,ip == 'BetaPF' &rif == 20 & rep == 1 & seq == "exo" &  isPeak == "yes")

df2[df2$summitPosDiff > 490,]

df2[df2$summitPosDiff < -500,]


df3 = subset(df,ip == 'Sig70' &rif == 0 & rep == 1 & seq == "exo" &  isPeak == "yes")

df3[df3$summitPosDiff > 500,]
df3[df3$summitPosDiff < -500,]

df4 = subset(df,ip == 'Sig70' &rif == 20 & rep == 1 & seq == "exo" &  isPeak == "yes")

df4[df4$summitPosDiff >400,]
df4[df4$summitPosDiff < -500,]
