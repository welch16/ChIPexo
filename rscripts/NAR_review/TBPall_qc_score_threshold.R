rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("TBP",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("chipseq",files,invert = TRUE)]

library(parallel)

files = files[-4]

reads = lapply(files,readGAlignments,param = NULL)
reads = lapply(reads,as,"GRanges")

set.seed(12345)

samp =c( seq_len(4) * 5e6 , seq(3,5) * 1e7)
lvs = c("5,000,000","10,000,000","15,000,000","20,000,000",
                                          "30,000,000","40,000,000","50,000,000")


samp_reads = lapply(reads,function(x){
    lapply(samp,function(y,z)z[sample(y)],x)})

names(samp_reads[[1]]) = prettyNum(as.integer(samp),big.mark = ",") 
names(samp_reads[[2]]) = prettyNum(as.integer(samp),big.mark = ",")
names(samp_reads[[3]]) = prettyNum(as.integer(samp),big.mark = ",")
names(samp_reads[[4]]) = prettyNum(as.integer(samp),big.mark = ",")

names(samp_reads) = c("Exo1","Exo2","Exo3","Nexus2")

samp_reads = unlist(samp_reads)

## small ntimes and nregions to be able to recalculate by ourselves
exo = lapply(samp_reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 22))

baseexo = lapply(reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 22))

library(data.table)

ntimes = 1e3
nregions = 1e3

dt = lapply(exo,as.data.frame)
dt = lapply(dt,as.data.table)

dt = lapply(dt,function(x)x[,.(uniquePos,depth)])

calculate_UParam <- function(stats,nregions,ntimes)
{

    calculate_UParam1 <- function(i,stats,nregions)
    {
        dt <- stats[sample(.N,nregions)]
        model <- lm(depth ~ 0 + uniquePos , data = dt)
        data.table(broom::tidy(model))
    }

    ss = mclapply(seq_len(ntimes),calculate_UParam1,stats,nregions,mc.cores = 20)
    rbindlist(ss)
   
}


library(dplyr)
library(tidyr)

totalReg = sapply(baseexo,length)

tt = tibble(exp = names(exo),nreg = sapply(exo,length))
tt = tt %>% separate(exp,into = c("repl","sample"),sep = "\\.")  %>%
    mutate(samp =as.numeric( gsub(",","",sample)),
           total = ifelse(repl == "Exo1",totalReg[1],
                   ifelse(repl == "Exo2",totalReg[2],
                   ifelse(repl == "Exo3",totalReg[3],totalReg[4]))),                   
           prop = nreg / total )

library(ggplot2)

figs = "figs/NAR_review/threshold"

theme_set(theme_bw())

pdf(file.path(figs,"TBPall_proportion_regions.pdf"))
u = tt %>% ggplot(aes(samp,prop,colour = repl))+geom_point()+geom_line()+ylim(0,1)
print(u)
dev.off()


beta1 = lapply(dt,calculate_UParam,nregions,ntimes)

beta1 = mapply(function(x,y)x[,name := y],beta1,names(beta1),SIMPLIFY = FALSE)


beta1 = rbindlist(beta1) %>% as.tbl %>% separate(name,into = c("repl","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = lvs))



pdf(file.path(figs,"TBPall_subsample_beta1.pdf"))
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("TBP in K562")
print(u)
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("TBP in K562")+scale_y_log10()
print(u)
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("TBP in K562")+ylim(0,800)
print(u)
dev.off()


bigdt = mapply(function(x,y)x[,name := y],dt,names(dt),SIMPLIFY = FALSE)
bigdt = rbindlist(bigdt) %>% as.tbl %>% separate(name,into = c("Rep","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = lvs),
           ratio = depth /uniquePos)


summaryRatio = bigdt %>% group_by(Rep,samp) %>%
    summarize(mean = mean(ratio),
              median = median(ratio),
              max = max(ratio),
              quant75 = quantile(ratio,prob = .75),
              quant90 = quantile(ratio,prob = .9),
              quant95 = quantile(ratio,prob = .95),
              quant99 = quantile(ratio,prob = .99))

summaryRatio %>% print
              

pdf(file.path(figs,"TBPall_depth_UniquePos_Ratio_histogram.pdf"),width = 12,height = 9)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-10,800)
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-15,40)
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-3,15)
print(u)
dev.off()
    



