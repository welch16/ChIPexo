
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

files = files[grep("vent",files)]

library(parallel)

reads = lapply(files,readGAlignments,param = NULL)
reads = lapply(reads,as,"GRanges")

set.seed(12345)

samp =c( seq_len(4) * 5e6 )#, seq(3,5) * 1e7)

samp_reads = lapply(reads,function(x){
    lapply(samp,function(y,z)z[sample(y)],x)})

names(samp_reads[[1]]) = prettyNum(as.integer(samp),big.mark = ",") 
names(samp_reads[[2]]) = prettyNum(as.integer(samp),big.mark = ",")
names(samp_reads[[3]]) = prettyNum(as.integer(samp),big.mark = ",") 

names(samp_reads) = c("Rep1","Rep2","Rep3")

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
           total = ifelse(repl == "Rep1",totalReg[1],
                   ifelse(repl == "Rep2",totalReg[2],totalReg[3])),
           prop = nreg / total
           ) 

library(ggplot2)

figs = "figs/NAR_review/threshold"

theme_set(theme_bw())

pdf(file.path(figs,"TBPexo_proportion_regions.pdf"))
tt %>% ggplot(aes(samp,prop,colour = repl))+geom_point()+geom_line()+ylim(0,1)
dev.off()


beta1 = lapply(dt,calculate_UParam,nregions,ntimes)

beta1 = mapply(function(x,y)x[,name := y],beta1,names(beta1),SIMPLIFY = FALSE)

beta1 = rbindlist(beta1) %>% as.tbl %>% separate(name,into = c("repl","samp"),sep = "\\.") %>%
    mutate(samp.nume = as.numeric(gsub(",","",samp)),
           samp = factor(samp, levels = c("5,000,000","10,000,000","15,000,000","20,000,000")))



pdf(file.path(figs,"TBPexo_subsample_beta1.pdf"))
beta1 %>% ggplot(aes(samp,estimate,fill = repl))+geom_boxplot()+facet_grid( . ~ repl)+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90))+
    geom_abline(slope = 0,intercept = c(0,10) ,linetype = 2)+
    ylab(expression(beta[1]))+ggtitle("ChIP-exo TBP in K562 ")+scale_y_log10()
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


##      Rep       samp      mean median       max quant75 quant90 quant95 quant99
##    <chr>     <fctr>     <dbl>  <dbl>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
## 1   Rep1  5,000,000  12.33975     12  247.5472  19.000    24.0      27   32.00
## 2   Rep1 10,000,000  12.42861     12  366.0000  19.000    24.0      27   32.00
## 3   Rep1 15,000,000  12.46539     12  366.0000  19.000    24.0      27   32.00
## 4   Rep1 20,000,000  12.51230     12  366.0000  19.000    24.0      27   32.00
## 5   Rep2  5,000,000 159.45366    111  859.0000 259.000   424.0     488  571.00
## 6   Rep2 10,000,000 157.91853    109  859.0000 254.000   422.0     488  571.00
## 7   Rep2 15,000,000 155.86559    105  964.0000 251.000   421.0     487  572.00
## 8   Rep2 20,000,000 155.48945    104  964.0000 251.000   422.0     488  572.00
## 9   Rep3  5,000,000 199.68166    150  988.0000 317.375   508.0     582  679.00
## 10  Rep3 10,000,000 195.58097    145 1241.3333 310.000   502.0     580  676.00
## 11  Rep3 15,000,000 194.31786    142 1241.3333 308.000   502.6     580  677.56
## 12  Rep3 20,000,000 193.06437    139 1241.3333 306.000   502.0     580  679.00


pdf(file.path(figs,"TBPexo_depth_UniquePos_Ratio_histogram.pdf"))
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 50,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")+
    xlim(-10,800)
print(u)
u = bigdt %>%
    ggplot(aes(ratio,fill = Rep))+
    geom_histogram(bins = 50,colour = "black",
                   aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
                   facet_grid(Rep ~ samp,scales = "free_y")+
    scale_fill_brewer(palette = "Pastel1",name = "Replicate")+
    xlab("D/U (0 to 40)")+theme(legend.position = "none")+ylab("Fraction of regions (normalized by panel)")+
    xlim(-4,40)
print(u)
dev.off()

                                     


bigdt

library(MASS)
