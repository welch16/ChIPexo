#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--outdr", action = "store_true",type = "character",
                default = "figs/figuresV2/fig5",
                help = "Directory where all the figures are stored."),
    make_option("--fig.width",action = "store_true",type = "numeric",default = 7,
                help = "Width of the figure composed with all panels"),
    make_option("--line.width",action = "store_true",type = "numeric",default = .8,
                help = "Plot line.width")    
)

opt = parse_args(OptionParser(option_list = optList))

library(grid)
library(gridExtra)
library(tidyverse)
library(magrittr)
library(viridis)
library(scales)

indr = "data/figures/fig5"



files = list.files(indr,full.names = TRUE,pattern = "tsv")

files = files[grep("scores",files)]

files1 = files[grep("M.tsv",files,invert = TRUE)]
files2 = files[grep("M.tsv",files,invert = FALSE)]

scores = lapply(files1,read_tsv)
scores = mapply(function(x,y)x %>% mutate(samp = y),scores,basename(files1),SIMPLIFY = FALSE) %>%
    bind_rows

scores = scores %>% mutate(samp = gsub("_scores.tsv","",samp)) %>%
    tidyr::separate(samp, into = c("lab","TF","Cell","Repl"),sep = "\\_") %>%
    mutate(protocol = ifelse(lab == "ChIPnexus","ChIP-nexus","ChIP-exo"),
           Cell = ifelse(Cell == "mouseliver","Mouse Liver",Cell),
           Cell = ifelse(Cell == "embryo","Embryo",Cell),
           Cell = ifelse(Cell == "yeast","Yeast",Cell),
           TF = ifelse(TF == "histone1","H3 Tail Deletion",TF)
           ) 



theme_set(theme_bw())

beta1 = scores %>%
    ggplot(aes(Repl,beta1,colour = protocol))+
    geom_boxplot(outlier.size = NA)+
    facet_grid( . ~ protocol + TF + Cell, space = "free_x",scales = "free_x")+ggtitle("A")+
    theme(legend.position = "top",
          axis.text.x = element_text(angle =90,vjust =.5),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 5))+
    scale_color_brewer(palette = "Dark2",name = "Protocol")+
    scale_y_log10(limits = c(1,400),breaks = c(1,10,50,100,500))+
    ylab(expression(beta[1]))

ggsave(file.path(opt$outdr,"fig5A.png"),width = 12,beta1)

beta2 = scores %>%
    ggplot(aes(Repl,-beta2,colour = protocol))+
    geom_boxplot(outlier.size = NA)+
    facet_grid( . ~ protocol + TF + Cell, space = "free_x",scales = "free_x")+ggtitle("B")+
    theme(legend.position = "none",
          axis.text.x = element_text(angle =90,vjust =.5),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 5))+
    scale_color_brewer(palette = "Dark2")+
    ylab(expression(beta[2]))+ylim(-.5,10)

ggsave(file.path(opt$outdr,"fig5B.png"),width = 12,beta2)






beta1_summary = scores %>%
    group_by(TF,Cell,lab,Repl) %>%
    summarize(
        N = n(),
        mean = mean(beta1),
        var = var(beta1),        
        median = median(beta1),        
        q1 = quantile(beta1,prob = .25),
        q3 = quantile(beta1,prob = .75),
        min = min(beta1),
        max = max(beta1),
        cor = cor(beta1,beta2)
    ) %>%
    mutate(IQR = q3 - q1,
           range = max - min)


beta2_summary = scores %>%
    group_by(TF,Cell,lab,Repl) %>%
    summarize(
        N = n(),
        mean = mean(beta2),
        var = var(beta2),        
        median = median(beta2),        
        q1 = quantile(beta2,prob = .25),
        q3 = quantile(beta2,prob = .75),
        min = min(beta2),
        max = max(beta2)                       
    ) %>%
    mutate(IQR = q3 - q1,
           range = max - min)


beta1_barplot = beta1_summary %>% ggplot(aes(Repl,log2(1 + IQR)))+geom_bar(stat = "identity")+
    facet_grid(paste(ifelse(grepl("Nexus",lab),"ChIP-nexus","ChIP-exo"),
                     TF, Cell)  ~ .,scales = "free",
               space = "free")+coord_flip()+
    theme(strip.text.y = element_text(angle = 0),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title.y = element_blank())+
    ylab(expression(log2(Delta ~ beta[1])))


beta2_barplot = 
beta2_summary %>% ggplot(aes(Repl,IQR))+geom_bar(stat = "identity")+
    facet_grid(paste(ifelse(grepl("Nexus",lab),"ChIP-nexus","ChIP-exo"),
                     TF, Cell)  ~ .,scales = "free",
               space = "free")+coord_flip()+
    theme(strip.text.y = element_text(angle = 0),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.title.y = element_blank())+
ylab(expression(Delta ~ beta[2]))



samp = lapply(files2,read_tsv)
samp = mapply(function(x,y)x %>% mutate(samp = y),samp,
              basename(files2),SIMPLIFY= FALSE) %>%
    bind_rows

samp = samp %>% mutate(samp = gsub("_scores","",samp),
                       samp = gsub(".tsv","",samp)) %>%
    separate(samp,into = c("lab","TF","Cell","repl","subsamp")) %>%
    mutate(prot = ifelse(lab == "venters","ChIP-exo","ChIP-nexus"),
           subsamp = as.numeric(gsub("M","",subsamp))) %>%
    mutate(TF = paste(prot , TF),
           Cell = ifelse(Cell == "mouseliver","Mouse Liver",Cell),
           Cell = ifelse(Cell == "embryo","Embryo",Cell),
           Cell = ifelse(Cell == "yeast","Yeast",Cell),
           TF = ifelse(TF == " histone1","H3 Tail Deletion",TF)
           ) 



samp_summary = samp %>% group_by(subsamp,prot,TF,repl) %>%
    summarize(beta1 = mean(beta1),
              beta2 = mean(beta2)
              )


beta1samp = samp_summary %>%
    ggplot(aes(subsamp,beta1,colour= prot,linetype = repl,shape = repl))+
    geom_line(show.legend = FALSE)+geom_point(size = 2)+
    scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
    scale_y_log10()+
    scale_color_brewer(palette = "Dark2", guide = FALSE)+
    ylab(expression(beta[1]))+
    xlab("Number of reads")+
    scale_shape_manual( values = c(0,1,2),name = "Replicate")+
    theme(legend.position = "top")+ggtitle("C")


g = ggplotGrob(beta1samp +theme(legend.position = "top"))$grobs
legend = g[[which(sapply(g,function(x)x$name) == "guide-box")]]

lh = sum(legend$height)


beta2samp = samp_summary %>%
    ggplot(aes(subsamp,-beta2,colour= prot,linetype = repl,shape = repl))+
    geom_line(show.legend = FALSE)+geom_point(size = 2)+
    scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
    scale_color_brewer(palette = "Dark2",
                       guide = FALSE)+
    theme(legend.position = "top")+
    ylab(expression(beta[2]))+
    scale_shape_manual( values = c(0,1,2),name = "Replicate")+
    ggtitle("D")+xlab("Number of reads")




ggsave(file.path(opt$outdr,"fig5C.png"),beta1samp)
ggsave(file.path(opt$outdr,"fig5D.png"),beta2samp)



fig5B = arrangeGrob(beta1samp + theme(legend.position = "none"),
                    beta2samp + theme(legend.position = "none"),nrow = 1)
fig5B = arrangeGrob(legend,fig5B,heights = unit.c(lh,unit(1,"npc") - lh),ncol =1)

fig5 = arrangeGrob(beta1,beta2,nrow = 2,heights = c(1.5,1))



all = arrangeGrob(beta1,beta2,fig5B,nrow = 3, heights = c(1.5,1,1))


ggsave(file = file.path(opt$outdr,"fig5.pdf"),all,width = 180,height = 220,units = "mm")






 




