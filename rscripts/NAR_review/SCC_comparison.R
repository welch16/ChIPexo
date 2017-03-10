
rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(magrittr)
library(ggplot2)

spp = "data/SCC_spp_curves"
chipqc = "data/SCC_curves"

files1 = list.files(spp,full.names = TRUE)
files2 = list.files(chipqc,full.names = TRUE)


files2 = sapply(files1,function(x)files2[grep(basename(x),files2)])
names(files2) = NULL

figs = "figs/NAR_review/SCC_comp"

scc_plot <- function(file1,file2)
{
    spp = read_delim(file1,delim = "\t")
    our = read_delim(file2,delim = "\t")
    scc = rbind(spp %>% mutate(program = "SPP"),
                our %>% mutate(program = "OUR")) %>%
        mutate(shift = as.numeric(shift))
    v = scc %>% group_by(program) %>% summarize(mm = which.max(cross.corr))
    scc %>% ggplot(aes(shift,cross.corr,colour = program)) + geom_line() +
        theme_bw()+ theme(legend.position = "top")+ scale_color_brewer(palette = "Set1") +
        ggtitle(gsub("_SCC.txt","",basename(file1)))+
        geom_vline(data = v,linetype = 2 ,aes(xintercept = mm,colour = program),show_guide = FALSE)
        
}

plots = mapply(scc_plot,files1,files2,SIMPLIFY = FALSE)


pdf(file.path(figs,"all_SCC_comparison.pdf"))
u = lapply(plots,print)
dev.off()
