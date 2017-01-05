
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)

dr = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(dr,full.names = TRUE)

files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]

library(parallel)

reads = mclapply(files,readGAlignments,param = NULL,mc.cores = 3)
reads = mclapply(reads,as,"GRanges",mc.cores = 3)

set.seed(12345)
samp = seq_len(4) * 5e6

samp_reads = mclapply(reads,function(x){
    lapply(samp,function(y,z)z[sample(y)],x)},mc.cores = 3)

names(samp_reads[[1]]) = prettyNum(as.integer(samp),big.mark = ",") 
names(samp_reads[[2]]) = prettyNum(as.integer(samp),big.mark = ",")
names(samp_reads[[3]]) = prettyNum(as.integer(samp),big.mark = ",") 

names(samp_reads) = c("Rep3","Rep1","Rep2")


samp_reads = unlist(samp_reads)

## small ntimes and nregions to be able to recalculate by ourselves
exo = lapply(samp_reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 20))

baseexo = lapply(reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 20))

## source("~/Desktop/Docs/Code/ChIPexoQual/R/base_summaryStats.R")

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
           total = ifelse(repl == "Rep1",totalReg[2],
                   ifelse(repl == "Rep2",totalReg[3],totalReg[1])),
           prop = nreg / total
           ) 

library(ggplot2)

figs = "figs/NAR_review/threshold"

theme_set(theme_bw())

pdf(file.path(figs,"FoxA1_proportion_regions.pdf"))
u = tt %>% ggplot(aes(samp,prop,colour = repl))+geom_point()+geom_line()+ylim(0,1)
print(u)
dev.off()


beta1 = lapply(dt,calculate_UParam,nregions,ntimes)

beta1 = mapply(function(x,y)x[,name := y],beta1,names(beta1),SIMPLIFY = FALSE)

beta1 = rbindlist(beta1) %>% as.tbl %>% separate(name,into = c("repl","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = c("5,000,000","10,000,000","15,000,000","20,000,000")))



pdf(file.path(figs,"FoxA1_subsample_beta1.pdf"))
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+ylim(0,12)+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("FoxA1 in mouse liver")
print(u)
dev.off()


bigdt = mapply(function(x,y)x[,name := y],dt,names(dt),SIMPLIFY = FALSE)
bigdt = rbindlist(bigdt) %>% as.tbl %>% separate(name,into = c("Rep","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = c("5,000,000","10,000,000","15,000,000","20,000,000")),
           ratio = depth /uniquePos)


summaryRatio = bigdt %>% group_by(Rep,samp) %>%
    summarize(mean = mean(ratio),
              median = median(ratio),
              max = max(ratio),
              quant75 = quantile(ratio,prob = .75),
              quant90 = quantile(ratio,prob = .9),
              quant95 = quantile(ratio,prob = .95),
              quant99 = quantile(ratio,prob = .99))

##      Rep       samp     mean median       max quant75 quant90 quant95 quant99
##    <chr>     <fctr>    <dbl>  <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
## 1   Rep1  5,000,000 1.443788      1  44.43038       2       2       3       4
## 2   Rep1 10,000,000 1.443261      1  44.43038       2       2       3       4
## 3   Rep1 15,000,000 1.442590      1 462.80208       2       2       3       4
## 4   Rep1 20,000,000 1.442829      1 462.80208       2       2       3       4
## 5   Rep2  5,000,000 1.222432      1  50.09783       1       2       2       3
## 6   Rep2 10,000,000 1.221634      1  50.09783       1       2       2       3
## 7   Rep2 15,000,000 1.221522      1 559.54497       1       2       2       3
## 8   Rep2 20,000,000 1.221540      1 559.54497       1       2       2       3
## 9   Rep3  5,000,000 6.555893      6  87.70000       9      12      14      17
## 10  Rep3 10,000,000 6.547870      6  87.70000       9      12      14      17
## 11  Rep3 15,000,000 6.552265      6 941.04790       9      12      14      17
## 12  Rep3 20,000,000 6.551857      6 941.04790       9      12      14      17



pdf(file.path(figs,"FoxA1_depth_UniquePos_Ratio_histogram.pdf"))
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    xlim(-2,25)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 10,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    xlim(-2,10)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
print(u)
dev.off()


## pdf(file.path(figs,"FoxA1_depth_UniquePos_Ratio_histogram1.pdf"))
## u = bigdt %>%
##     filter(depth > 10) %>%
##     ggplot(aes(ratio,fill = Rep))+
##     geom_histogram(bins = 25,colour = "black",
##                    aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
##                    facet_grid(Rep ~ samp,scales = "free_y")+
##     xlim(-2,25)+
##     scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
##     xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
## print(u)
## u = bigdt %>%
##     filter(depth > 10) %>%    
##     ggplot(aes(ratio,fill = Rep))+
##     geom_histogram(bins = 10,colour = "black",
##                    aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
##                    facet_grid(Rep ~ samp,scales = "free_y")+
##     xlim(-2,10)+
##     scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
##     xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
## print(u)
## dev.off()




ex = bigdt %>% filter(Rep == "Rep3") %>% group_by(samp)


