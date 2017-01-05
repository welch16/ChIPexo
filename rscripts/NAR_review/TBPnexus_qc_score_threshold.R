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
files = files[grep("vent",files,invert = TRUE)]

library(parallel)

reads = lapply(files,readGAlignments,param = NULL)
reads = lapply(reads,as,"GRanges")

set.seed(12345)

samp =c( seq_len(6) * 5e6)

samp_reads = lapply(reads,function(x){
    lapply(samp,function(y,z)z[sample(y)],x)})

names(samp_reads[[1]]) = prettyNum(as.integer(samp),big.mark = ",") 
names(samp_reads[[2]]) = prettyNum(as.integer(samp),big.mark = ",")

## names(samp_reads[[3]]) = prettyNum(as.integer(samp),big.mark = ",") 

names(samp_reads) = c("Rep1","Rep2")

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
           total = ifelse(repl == "Rep1",totalReg[1],totalReg[2]),
           prop = nreg / total
           )

library(ggplot2)

figs = "figs/NAR_review/threshold"

theme_set(theme_bw())

pdf(file.path(figs,"TBPnexus_proportion_regions.pdf"))
tt %>% ggplot(aes(samp,prop,colour = repl))+geom_point()+geom_line()+ylim(0,1)
dev.off()


beta1 = lapply(dt,calculate_UParam,nregions,ntimes)

beta1 = mapply(function(x,y)x[,name := y],beta1,names(beta1),SIMPLIFY = FALSE)

lvs =  c("5,000,000","10,000,000","15,000,000","20,000,000","25,000,000","30,000,000")

beta1 = rbindlist(beta1) %>% as.tbl %>% separate(name,into = c("repl","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = lvs))



pdf(file.path(figs,"TBPnexus_subsample_beta1.pdf"))
beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+ylim(0,40)+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("ChIPnexus TBP in K562")
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


##      Rep       samp      mean median      max quant75 quant90 quant95 quant99
##    <chr>     <fctr>     <dbl>  <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
## 1   Rep1  5,000,000  2.557556      2 211.5960     3.0       5       6       8
## 2   Rep1 10,000,000  2.532458      2 211.5960     3.0       5       6       8
## 3   Rep1 15,000,000  2.530728      2 250.8311     3.0       5       6       8
## 4   Rep1 20,000,000  2.533460      2 250.8311     3.0       5       6       8
## 5   Rep1 25,000,000  2.535610      2 250.8311     3.0       5       6       8
## 6   Rep1 30,000,000  2.543696      2 250.8311     3.0       5       6       8
## 7   Rep2  5,000,000 16.031849     15 877.1515    21.0      30      34      41
## 8   Rep2 10,000,000 15.602061     14 877.1515    20.5      29      34      41
## 9   Rep2 15,000,000 15.474208     14 877.1515    20.0      29      34      41
## 10  Rep2 20,000,000 15.441075     14 877.1515    20.0      29      34      41
## 11  Rep2 25,000,000 15.342682     14 877.1515    20.0      29      33      41
## 12  Rep2 30,000,000 15.305654     14 877.1515    20.0      29      33      41


pdf(file.path(figs,"TBPnexus_depth_UniquePos_Ratio_histogram.pdf"),width = 10)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 50,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+xlim(-5,50)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 15,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+xlim(-2,15)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")
print(u)
dev.off()
    
