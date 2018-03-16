#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--outdr", action = "store_true",type = "character",
                default = "figs/figuresV2/fig1/thesis",
                help = "Directory where all the figures are stored."),
    make_option("--fig.width",action = "store_true",type = "numeric",default = 7,
                help = "Width of the figure composed with all panels"),
    make_option("--line.width",action = "store_true",type = "numeric",default = .8,
                help = "Plot line.width")    
)

opt = parse_args(OptionParser(option_list = optList))


library(gridExtra)
library(tidyverse)
library(magrittr)
library(viridis)
library(scales)

indr = "data/figures/fig1"

files = list.files(indr,full.names = TRUE,pattern = "tsv")


theme_set(theme_bw())

r = viridis(100,option = "D")

bindata = files[grep("bins",files)] %>% read_tsv

chipexo_chipseqPE_comparison = bindata %>%
    ggplot(aes(ChIPseqPET,ChIPexo))+stat_binhex(bins = 100)+
    xlab("ChIP-seq (PE) counts")+ylab("ChIP-exo counts")+
    scale_x_log10()+scale_y_log10()+
    scale_fill_gradientn(colours = r,trans = 'log10',
                         labels=trans_format('log10',math_format(10^.x)),
                         name = "Number of bins")+
    theme(legend.position = "top")+
    ggtitle("A")

dir.create(opt$outdr,showWarnings = FALSE)

## mappability

mapdata = files[grep("map",files)] %>% read_tsv %>%
    mutate(meanSq = varYall + meanYall^2)

tiles = 50
z = 1.96
mapdata2 = mapdata %>% mutate(groups = ntile(uS,tiles) / tiles) %>%
    group_by(groups,sample) %>%
    summarize(uS = median(uS),
              allnitem = sum(nitem),
              meanYall = weighted.mean(meanYall , w = nitem / allnitem ),
              meanSqall = weighted.mean(meanSq, w= nitem / allnitem,na.rm = TRUE),
              varYall = meanSqall - meanYall^2
              ) %>%
    mutate( lb = meanYall - z * sqrt(varYall / allnitem),
           ub = meanYall + z * sqrt(varYall / allnitem)) %>%
    filter(!(between(uS,.3,.5) & meanYall <.5))



mapplot = mapdata2 %>%
    ggplot(aes(uS,meanYall,color = sample))+geom_point()+
    geom_linerange(data = mapdata2,aes(group =1,ymin = lb,ymax = ub),
                   size = opt$line.width,show.legend = FALSE)+
    theme(legend.position = "top")+
    scale_color_brewer(palette = "Dark2",name = "Replicate")+
    xlab("Mappability Score")+ylab("Average ChIP read counts")+    
    ggtitle("B")

## GC content

gcdata = files[grep("GC",files)] %>% read_tsv %>%
    mutate(meanSq = varYall + meanYall^2)

gcdata2 = gcdata %>% mutate(groups = ntile(uS,tiles) / tiles) %>%
    group_by(groups,sample) %>%
    summarize(uS = median(uS),
              allnitem = sum(nitem),
              meanYall = weighted.mean(meanYall,w = nitem / allnitem),
              meanSqall = weighted.mean(meanSq,w = nitem / allnitem,na.rm = TRUE),
              varYall = meanSqall - meanYall^2) %>%
    mutate(
        lb = meanYall - z * sqrt(varYall / allnitem),
        ub = meanYall + z * sqrt(varYall / allnitem),
        range = ub - lb) %>% filter(range < 3) %>%
    filter(!( between(uS,.5,.6) & meanYall < 1))

gcplot = gcdata2 %>% ggplot(aes(uS,meanYall,color = sample))+
    geom_point()+
    geom_linerange(data = gcdata2,aes(ymin = lb,ymax = ub,group = 1),
                   size = opt$line.width,show.legend = FALSE)+
    theme(legend.position = "top")+
    scale_color_brewer(palette = "Dark2",name = "Replicate")+
    xlab("GC content Score")+ylab("Average ChIP read counts")+ 
    ggtitle("C")


## SCC

scc = files[grep("scc",files)] %>% read_tsv


scc = scc %>% mutate(prot = ifelse(grepl("exo",sample),"ChIP-exo","ChIP-seq (SE)"),
                     repl = ifelse(prot == "ChIP-exo",gsub("ChIP-exo","",sample),
                                   gsub("ChIP-seq","",sample)),
                     repl = paste0("Rep",repl))

rr = brewer_pal(4,palette = "Set1")(4)[c(1,3)]

sccplot = scc %>% ggplot(aes(shift,cross.corr,colour = sample))+
    geom_line(size = opt$line.width)+
    facet_grid( . ~ prot  )+geom_vline(xintercept = 36,linetype = 3,colour = "grey")+
    theme(legend.position = "top")+ggtitle("D")+
    xlab("Shift")+ylab("Strand Cross-Correlation")+
    scale_color_brewer(palette = "Dark2",name = "",
                       guide = guide_legend(nrow = 2))
    

all = arrangeGrob(
             chipexo_chipseqPE_comparison,
             mapplot,
             gcplot,
             sccplot,
             ncol = 2,nrow = 2,heights = c(1.5,1.5)
             )

ggsave(file = file.path(opt$outdr,"fig2_thesis.png"),all,width = 180,height = 180,units = "mm",
       dpi = 1400)



