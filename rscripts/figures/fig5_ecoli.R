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

indr = "data/figures/fig5/ecoli/800"



files = list.files(indr,full.names = TRUE,pattern = "tsv")

files = files[grep("scores",files)]

files1 = files[grep("M.tsv",files,invert = TRUE)]
files2 = files[grep("M.tsv",files,invert = FALSE)]


scores = lapply(files1,read_tsv)
scores = mapply(function(x,y)x %>% mutate(samp = y),
                scores,basename(files1),SIMPLIFY = FALSE) %>%
    bind_rows()

scores = scores %>%
    mutate(samp = gsub("_scores.tsv","",samp)) %>%
    tidyr::separate(samp, into = c("lab","sample"),sep = "\\_") %>%
    mutate(TF = "Sig70",
           Growth = ifelse(sample %in% c("935","937"),"Stationary","Exponential"),
           Rif = ifelse(sample %in% c("1314","1320"),"Rif_20min","Rif_0min"),
           rep = ifelse(sample %in% c("931","935","1311","1314"),"Rep-1","Rep-2"),
           group = ifelse(sample %in% c("931","935","933","937"),"E1","E2"),
           Repl = ifelse(group == "E1",paste(Growth,rep,sep = "_"),
                         paste(Rif,rep,sep = "_"))
           )



theme_set(theme_bw())

beta1 = scores %>%
    ggplot(aes(Repl,beta1,colour = group))+
    geom_boxplot(outlier.size = NULL)+
    facet_grid( . ~ group, space = "free_x",scales = "free_x")+ggtitle("A")+
    theme(legend.position = "top",
          axis.text.x = element_text(angle =90,vjust =.5),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 5))+
    scale_color_brewer(palette = "Set1",name = "Group")+
    scale_y_log10(limits = c(1,400),breaks = c(1,10,50,100,400))+    
    ylab(expression(beta[1]))

ggsave(file.path(opt$outdr,"fig5A_ecoli.png"),width = 12,beta1)

beta2 = scores %>%
    ggplot(aes(Repl,-beta2,colour = group))+
    geom_boxplot(outlier.size = NULL)+
    facet_grid( . ~ group, space = "free_x",scales = "free_x")+ggtitle("B")+
    scale_color_brewer(palette = "Set1",name = "Group")+    
    theme(legend.position = "none",
          axis.text.x = element_text(angle =90,vjust =.5),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 5))+
    ylab(expression(beta[2]))

ggsave(file.path(opt$outdr,"fig5B_ecoli.png"),width = 12,beta2)






beta1_summary = scores %>%
    group_by(group,Repl) %>%
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
    group_by(group,Repl) %>%
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



samp = map(files2,read_tsv) %>%
    map2(basename(files2),
            function(x,y)x %>% mutate(samp = y)) %>%
    bind_rows()

samp = samp %>% mutate(samp = gsub("_scores_","",samp),
                       samp = gsub(".tsv","",samp)) %>%
    separate(samp , into = c("sample","TF","subsamp_lab"),sep = "\\_") %>%
    mutate(TF = "Sig70",
           sample = gsub("edsen","",sample),
           Growth = ifelse(sample %in% c("935","937"),"Stationary","Exponential"),
           Rif = ifelse(sample %in% c("1314","1320"),"Rif_20min","Rif_0min"),
           rep = ifelse(sample %in% c("931","935","1311","1314"),"Rep-1","Rep-2"),
           group = ifelse(sample %in% c("931","935","933","937"),"E1","E2"),
           Repl = ifelse(group == "E1",paste(Growth,rep,sep = "_"),
                         paste(Rif,rep,sep = "_")))%>%
    mutate(subsamp = gsub("subsample","",subsamp_lab),
           subsamp = 1e6 * as.numeric(gsub("M","",subsamp)),
           subsamp_lab = paste0(subsamp / 1e3,"K"))%>%
    mutate(
        Repl = factor(Repl,levels = c("Exponential_Rep-1",
                                      "Exponential_Rep-2",
                                      "Stationary_Rep-1",
                                      "Stationary_Rep-2",
                                      "Rif_0min_Rep-1",
                                      "Rif_0min_Rep-2",
                                      "Rif_20min_Rep-1",
                                      "Rif_20min_Rep-2")))
                                      
                                      
    
           
samp_summary = samp %>% group_by(group,subsamp,Repl) %>%
    summarize(beta1 = mean(beta1),
              beta2 = mean(beta2)
              )


beta1samp = samp_summary %>%
    ggplot(aes(subsamp,beta1,colour =group,linetype = Repl,shape = Repl))+
    geom_line(show.legend = FALSE)+geom_point(size = 2)+
    scale_color_brewer(palette = "Set1", guide = FALSE)+
    ylab(expression(beta[1]))+
    xlab("Number of reads")+
    scale_shape_manual(name = "Sample",values = rep(c(0,1,2,4),2))+
    theme(legend.position = "top")+ggtitle("C")


g = ggplotGrob(beta1samp +theme(legend.position = "top"))$grobs
legend = g[[which(sapply(g,function(x)x$name) == "guide-box")]]

lh = sum(legend$height)

beta2samp = samp_summary %>%
    ggplot(aes(subsamp,-beta2,colour =group,linetype = Repl,shape = Repl))+
    geom_line(show.legend = FALSE)+geom_point(size = 2)+
    scale_color_brewer(palette = "Set1", guide = FALSE)+
    ylab(expression(beta[2]))+
    xlab("Number of reads")+
    scale_shape_manual(name = "Sample",values = rep(c(0,1,2,4),2))+
    theme(legend.position = "top")+ggtitle("D")




ggsave(file.path(opt$outdr,"fig5C_ecoli.png"),beta1samp)
ggsave(file.path(opt$outdr,"fig5D_ecoli.png"),beta2samp)



fig5B = arrangeGrob(beta1samp + theme(legend.position = "none"),
                    beta2samp + theme(legend.position = "none"),nrow = 1)
fig5B = arrangeGrob(legend,fig5B,heights = unit.c(lh,unit(1,"npc") - lh),ncol =1)

fig5 = arrangeGrob(beta1,beta2,nrow = 2,heights = c(1.5,1))



all = arrangeGrob(beta1,beta2,fig5B,nrow = 3, heights = c(1.5,1,1))


ggsave(file = file.path(opt$outdr,"fig5_ecoli.pdf"),all,width = 180,height = 300,units = "mm")







 




