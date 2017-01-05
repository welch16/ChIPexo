rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"
files = list.files(dr,full.names = TRUE,recursive = TRUE)

files = files[grep("mei",files)]
files = files[grep("sorlt",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("U2OS",files,invert = TRUE)]

library(parallel)

reads = lapply(files,readGAlignments,param = NULL)
reads = lapply(reads,as,"GRanges")

set.seed(12345)

samp =c( seq_len(4) * 5e6 , seq(3,4) * 1e7)
lvs = c("5,000,000","10,000,000","15,000,000","20,000,000",
                                          "30,000,000","40,000,000")

samp_reads = lapply(reads,function(x){
    lapply(samp,function(y,z)z[sample(y)],x)})

names(samp_reads[[1]]) = prettyNum(as.integer(samp),big.mark = ",") 
names(samp_reads[[2]]) = prettyNum(as.integer(samp),big.mark = ",")

names(samp_reads) = c("IMR90","K562")

samp_reads = unlist(samp_reads)

## small ntimes and nregions to be able to recalculate by ourselves
exo = lapply(samp_reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 22))

## baseexo = lapply(reads,function(x)ExoData(reads = x, nregions = 1000,ntimes = 100,mc.cores = 22))

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

## totalReg = sapply(baseexo,length)

## tt = tibble(exp = names(exo),nreg = sapply(exo,length))
## tt = tt %>% separate(exp,into = c("repl","sample"),sep = "\\.")  %>%
##     mutate(samp =as.numeric( gsub(",","",sample)),
##            total = ifelse(repl == "Exo1",totalReg[1],
##                    ifelse(repl == "Exo2",totalReg[2],
##                    ifelse(repl == "Exo3",totalReg[3],totalReg[4]))),                   
##            prop = nreg / total )

library(ggplot2)

figs = "figs/NAR_review/threshold"

theme_set(theme_bw())

## pdf(file.path(figs,"TBPall_proportion_regions.pdf"))
## u = tt %>% ggplot(aes(samp,prop,colour = repl))+geom_point()+geom_line()+ylim(0,1)
## print(u)
## dev.off()


beta1 = lapply(dt,calculate_UParam,nregions,ntimes)

beta1 = mapply(function(x,y)x[,name := y],beta1,names(beta1),SIMPLIFY = FALSE)


beta1 = rbindlist(beta1) %>% as.tbl %>% separate(name,into = c("repl","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = lvs))



pdf(file.path(figs,"GR_subsample_beta1.pdf"))
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("GR")
print(u)
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("GR")+scale_y_log10()
print(u)
u = beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("GR")+ylim(0,30)
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

##      Rep       samp      mean   median       max quant75 quant90 quant95
##    <chr>     <fctr>     <dbl>    <dbl>     <dbl>   <dbl>   <dbl>   <dbl>
## 1  IMR90  5,000,000  6.596513  5.00000  66.00000    10.0    15.0      19
## 2  IMR90 10,000,000  6.615544  5.00000  66.00000    10.0    15.0      19
## 3  IMR90 15,000,000  6.591771  5.00000  66.00000    10.0    15.0      19
## 4  IMR90 20,000,000  6.591833  5.00000  66.00000    10.0    15.0      19
## 5  IMR90 30,000,000  6.541999  5.00000  73.83750    10.0    15.0      19
## 6  IMR90 40,000,000  6.544241  5.00000  73.83750    10.0    15.0      19
## 7   K562  5,000,000 14.323767 13.00000  83.36842    19.0    26.5      32
## 8   K562 10,000,000 13.940911 13.00000 313.86463    19.0    26.0      32
## 9   K562 15,000,000 13.886006 12.66667 313.86463    19.0    26.0      32
## 10  K562 20,000,000 13.793588 12.40000 313.86463    18.5    26.0      32
## 11  K562 30,000,000 13.710742 12.00000 313.86463    18.5    26.0      32
## 12  K562 40,000,000 13.577788 12.00000 313.86463    18.0    26.0      32


pdf(file.path(figs,"GR_depth_UniquePos_Ratio_histogram.pdf"),width = 12,height = 9)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-3,30)
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-10,50)
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 25,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions")+xlim(-10,100)
print(u)
dev.off()
    



