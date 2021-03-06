#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--outdr", action = "store_true",type = "character",
                default = "figs/figuresV2/fig2",
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

indr = "data/figures/fig3"

files = list.files(indr,full.names = TRUE,pattern = "tsv")


theme_set(theme_bw())

r = viridis(100,option = "D")

dir.create(opt$outdr,showWarnings = FALSE)

##
exo = files[grep("ARCv",files)] %>% read_tsv

exo = exo %>% filter(Repl == "Rep-1")


ARC = exo %>% ggplot(aes(ARC,URC))+stat_binhex(bins = 50)+
    scale_fill_gradientn(colours = r , trans = "log10",
                         labels = trans_format("log10",
                                               math_format(10^.x)),
                         name = "Nr. of islands")+
    theme(legend.position = "top") + facet_grid( ~ Repl)+xlim(0,4)


ggsave(file.path(opt$outdr,"fig2A.png"),ARC)

## Region Comp

comp = files[grep("RegionComp",files)] %>% read_tsv

comp = comp %>% filter(Repl == "Rep-1")

r = brewer_pal(palette = "Pastel1")

region_comp = comp %>%
    ggplot(aes(depth,prob,fill = lab))+geom_bar(stat = "identity")+
    facet_grid(  Repl ~ .)+
    theme(legend.position = "top")+scale_fill_manual(values = r(3),name = "")+
    ylab("Proportion of islands")+
    xlab("Min. number of reads")
        

ggsave(file.path(opt$outdr,"fig2B.png"),width = 6 ,height = 4,region_comp)

## FSR dist

FSRdist = files[grep("FSR",files)] %>% read_tsv

FSRdist = FSRdist %>% filter(Repl == "Rep-1")


FSR = FSRdist %>% ggplot(aes(depth,FSR,colour =as.factor( quantiles)))+
    geom_line(,size =opt$line.width)+
    theme(legend.position = "top")+facet_grid( Repl ~ .)+
    ylab("Forward Strand Ratio (FSR)")+xlab("Min. number of reads")+
    scale_color_brewer(palette = "Set1",name = "Quantiles",
                       guide = guide_legend(nrow = 2))


ggsave(file.path(opt$outdr,"fig2C.png"),width = 6,height = 4 , FSR)


## Blacklist boxplot

blacklists = files[grep("blacklist",files)] %>% read_tsv %>%
    mutate(repl = forcats::fct_recode(repl,"Rep-3"="rep-3","Rep-2"="rep-2","Rep-1"="rep-1"))

library(ggthemes)


blacklists = blacklists %>% filter(repl == "Rep-1")

beta1 = blacklists %>% filter(term == "uniquePos") %>%    
    ggplot(aes(repl,estimate,colour = blackl))+geom_boxplot()+
    theme(legend.position = "top",
          axis.title.x = element_blank())+
    scale_color_brewer(palette = "Set1",name = "")+
    scale_y_log10()+ylab(expression(beta[1]))

beta2 = blacklists %>% filter(term != "uniquePos") %>%    
    ggplot(aes(repl,-estimate,colour = blackl))+geom_boxplot()+
    scale_color_brewer(palette = "Set1",name = "")+    
    theme(legend.position = "top",
          axis.title.x = element_blank())+
    ylab(expression(beta[2]))+ylim(0,10)

grid_arrange_shared_legend <- function(...)
{
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = "top"))$grobs
    legend <- g[[which(sapply(g, function(x)x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        legend,
        do.call(arrangeGrob,lapply(plots,function(x)
            x + theme(legend.position = "none"))),        
        ncol = 1,
        heights = unit.c(unit(1,"npc") - lheight,lheight))   
}


g = ggplotGrob(beta1 +theme(legend.position = "top"))$grobs
legend = g[[which(sapply(g,function(x)x$name) == "guide-box")]]

lh = sum(legend$height)

betascore = arrangeGrob(beta1 + theme(legend.position = "none"),
                        beta2 + theme(legend.position = "none"),nrow = 1)

betascore = arrangeGrob(legend,betascore,heights = unit.c(lh,unit(1,"npc") - lh),ncol =1)

ggsave(file.path(opt$outdr,"fig2D.png"),width = 6,betascore)


all = arrangeGrob(arrangeGrob(ARC,nrow = 1,ncol = 1),
                  arrangeGrob(region_comp,FSR,nrow = 1,ncol = 2),
                  betascore,
                  heights = c(1.2,1.3,.85))

ggsave(file = file.path(opt$outdr,"fig2.pdf"),all,
       width = 180,
       units = "mm")
