#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--outdr", action = "store_true",type = "character",
                default = "figs/figuresV2/fig4",
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

indr = "data/figures/fig4"

files = list.files(indr,full.names = TRUE,pattern = "tsv")


theme_set(theme_bw())

dir.create(opt$outdr,showWarnings = FALSE)

## FIMO SCORES  - FoxA1

FOXA1_scores = files[grepl("FoxA1",files) &
                     grepl("FIMO",files)] %>% read_tsv

foxA1_boxplot = FOXA1_scores %>%
    ggplot(aes(Replicate,score,fill = Replicate))+geom_boxplot()+
    facet_grid( . ~ K)+scale_fill_brewer(palette = "Pastel1")+
    ylab("FIMO scores")+
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))+
    ggtitle("A")

ggsave(file.path(opt$outdr,"fig4A.png"),foxA1_boxplot)

## Profiles

prof = files[grepl("profiles",files)] %>% read_tsv %>%
    mutate(which = paste0(strand,":",rep),
           strand = forcats::fct_recode(strand,
                                 "Forward" = "+",
                                 "Reverse" = "-"
                                 ))
                                

r4 = c("darkblue","firebrick3","blue","red","lightblue","lightpink")

profiles =  prof %>%
    ggplot(aes(coord,counts,colour = which))+geom_line()+
    scale_colour_manual(values = r4,
                        labels = c("Reverse Rep-1","Forward Rep-1","Reverse Rep-2",
                                   "Forward Rep-2","Reverse Rep-3","Forward Rep-3"),name = "")+
    theme(legend.position = "top",
          legend.text = element_text(size = 6))+
    xlab("Position around motif start")+
    ylab("Average counts")+ggtitle("B")


ggsave(file.path(opt$outdr,"fig4B.png"),profiles)

## FoxA1 peak pair assumption

FoxA1_peakPair = files[grepl("FoxA1",files) &
                       grepl("imbalance",files)] %>% read_tsv

FSR_hist = FoxA1_peakPair %>%
    ggplot(aes(FSR,fill = nmotif)) + geom_histogram(colour = "black",bins = 25)+
    facet_grid( ~ repl)+
    theme(legend.position = "top",
          legend.text = element_text(size = 6))+
    scale_fill_brewer(palette = "Pastel2",
                      name = "Nr. of motifs")+
    ylab("ChIP regions overlapping peaks")+
    xlab("Forward Strand Ratio (FSR)")+
    scale_x_continuous(breaks = c(0,.5,1))
    
ggsave(file.path(opt$outdr,"fig4C.png"),FSR_hist)    


## TBP fimo scores


TBP_scores = files[grepl("TBP",files) &
                   grepl("FIMO",files)] %>% read_tsv %>%
             mutate(aux = replicate,
                    replicate = gsub("Exo","ChIP-exo::",replicate),
                    replicate = gsub("Nexus","ChIP-nexus::",replicate)
                    ) %>%
             tidyr::separate(replicate,into = c("prot","repl"),sep ="::") %>%
             mutate(repl = paste0("Rep-",repl)) 
             

TBP_boxplot = TBP_scores %>%
    ggplot(aes(aux,score,fill = prot))+geom_boxplot()+
    facet_grid( . ~ K)+scale_fill_brewer(palette = "Pastel2",name = "")+
    ylab("FIMO scores")+
    theme(
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))+
    scale_x_discrete(labels = paste0("Rep-",c(1,2,3,1,2)))+
    ggtitle("D")

ggsave(file.path(opt$outdr,"fig4D.png"),TBP_boxplot)


## FoxA1 peak pair assumption

TBP_peakPair = files[grepl("TBP",files) &
                     grepl("imbalance",files)] %>% read_tsv %>%
               mutate(nmotif =
                          forcats::fct_collapse(as.factor(nmotif),">=3" = as.character(seq(3,100))),
                      repl = forcats::fct_relevel(repl,
                                                  "ChIP-exo-1"="Exo1",
                                                  "ChIP-exo-2"="Exo2",
                                                  "ChIP-exo-3"="Exo3",
                                                  "ChIP-nexus-1"="Nexus1",
                                                  "ChIP-nexus-2"="Nexus2"))


                                                  
FSR_hist21 = TBP_peakPair %>% filter(!grepl("Nexus",repl)) %>%
    ggplot(aes(FSR,fill = as.factor(nmotif))) + geom_histogram(colour = "black",bins = 25)+
    facet_grid( ~ repl)+
    theme(legend.position = "top",legend.text = element_text(size = 6))+
    scale_fill_brewer(palette = "Pastel2",
                      name = "Nr. of motis")+
    ylab("ChIP regions overlapping peaks")+
    xlab("Forward Strand Ratio (FSR)")+
    scale_x_continuous(breaks = c(0,.5,1))


FSR_hist22 = TBP_peakPair %>% filter(grepl("Nexus",repl)) %>%
    ggplot(aes(FSR,fill = as.factor(nmotif))) + geom_histogram(colour = "black",bins = 25)+
    facet_grid( ~ repl)+
    theme(legend.position = "top",legend.text = element_text(size = 6))+
    scale_fill_brewer(palette = "Pastel2",
                      name = "Nr. of motifs")+
    ylab("ChIP regions overlapping peaks")+
    xlab("Forward Strand Ratio (FSR)")+
    scale_x_continuous(breaks = c(0,.5,1))


g = ggplotGrob(FSR_hist21 +theme(legend.position = "top"))$grobs
legend = g[[which(sapply(g,function(x)x$name) == "guide-box")]]

lh = sum(legend$height)

TBP_FSR = arrangeGrob(FSR_hist21 + theme(legend.position = "none") +
                     ggtitle("E"),
                     FSR_hist22 + theme(legend.position = "none")+
                     ggtitle("") +
                    theme(axis.title.y = element_blank()),nrow = 1)

TBP_FSR = arrangeGrob(legend,TBP_FSR,heights = unit.c(lh,unit(1,"npc") - lh),ncol =1)


ggsave(file.path(opt$outdr,"fig4E.png"),TBP_FSR,width = 12)    


all = arrangeGrob(arrangeGrob(foxA1_boxplot,nrow = 1,ncol = 1),
                  arrangeGrob(profiles,FSR_hist,nrow = 1,ncol = 2),
                  arrangeGrob(TBP_boxplot),
                  TBP_FSR,
                  heights = c(.6,1,.8,.8))

ggsave(file = file.path(opt$outdr,"fig4.pdf"),all,
       width = 180,
       height = 280,
       units = "mm")

