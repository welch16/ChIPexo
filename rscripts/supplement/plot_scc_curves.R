
rm(list = ls())

library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)

data_dir <- "data/ChIPseq_QC"
files <- list.files(data_dir,recursive = TRUE,full.names = TRUE,
                    include.dirs = TRUE)

qc_ind <- lapply(files,read.table,header = TRUE)
qc_ind <- lapply(qc_ind,data.table)

qc_ind <- mapply(function(x,y)x[,name := y],qc_ind,
                 gsub("_QC.txt","",basename(files)),
                 SIMPLIFY = FALSE)
qc_ind <- do.call(rbind,qc_ind)
qc_ind <- qc_ind %>% filter(read_length < 100)

dr <- "data/SCC_curves"
figs_dir <- "figs/supplement"

files <- list.files(dr,full.names = TRUE,include.dirs = TRUE)
nms <- basename(files)

scc <- lapply(files,fread)

scc <- mapply(function(x,y)x[,name := gsub("_SCC.txt","",y)],
   scc,nms,SIMPLIFY = FALSE)
scc <- do.call(rbind,scc)


r1 <- brewer.pal(name = "Set1",5)
r2 <- brewer.pal(name = "Dark2",5)
r <- c(r2[c(1,3:5)],r1[2:3])

## carroll human ## need to compare to CHIP-Seq

rl <- qc_ind %>% filter(grepl("carroll_human",name)) %>%
  select(read_length) %>% unique()

carroll_human <- scc %>% filter(grepl("carroll_human",name) |
                                grepl("chipseq",name)) %>%
  mutate(cond = ifelse(grepl("chipseq",name),"SE ChIP-Seq","ChIP-exo"),
         repl = ifelse(grepl("rep1",tolower(name)),"Rep-1",
           ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Carroll_human_SCC.pdf"),width = 6, height = 4)
ggplot(carroll_human,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  facet_grid( . ~ cond)+
  geom_vline(xintercept = rl[[1]],linetype = 2)+
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()

## carroll mouse

rl <- qc_ind %>% filter(grepl("carroll_mouse",name)) %>%
  select(read_length) %>% unique()


carroll_mouse <- scc %>% filter(grepl("carroll_mouse",name)) %>%
  mutate(repl = ifelse(grepl("rep1",tolower(name)),"Rep-1",
           ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Carroll_mouse_SCC.pdf"),width = 4, height = 4)
ggplot(carroll_mouse,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+  
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()

## nexus embryo dorsal

rl <- qc_ind %>% filter(grepl("embryo_dorsal",name)) %>%
  select(read_length) %>% unique()

embryo_dorsal <- scc %>% filter(grepl("embryo_dorsal",name)) %>%
  mutate(repl = ifelse(grepl("rep1",tolower(name)),"Rep-1",
           ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Nexus_embryo_dorsal_SCC.pdf"),width = 4, height = 4)
ggplot(embryo_dorsal,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+  
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()

## nexus embryo twist
rl <- qc_ind %>% filter(grepl("embryo_twist",name)) %>%
  select(read_length) %>% unique()


embryo_twist <- scc %>% filter(grepl("embryo_twist",name)) %>%
  mutate(repl = ifelse(grepl("rep1",tolower(name)),"Rep-1",
           ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Nexus_embryo_twist_SCC.pdf"),width = 4, height = 4)
ggplot(embryo_twist,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()

## nexus s2 max
rl <- qc_ind %>% filter(grepl("S2_Max",name)) %>%
  select(read_length) %>% unique()

s2_max <- scc %>% filter(grepl("S2_Max",name)) %>% mutate(repl =
  ifelse(grepl("rep1",tolower(name)),"Rep-1",
  ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Nexus_S2_Max_SCC.pdf"),width = 4,height = 4)
ggplot(s2_max,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+ theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+  
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")
dev.off()

## nexus s2 myc
rl <- qc_ind %>% filter(grepl("S2_MyC",name)) %>%
  select(read_length) %>% unique()


s2_myc <- scc %>% filter(grepl("S2_MyC",name)) %>% mutate(repl =
  ifelse(grepl("rep1",tolower(name)),"Rep-1",
  ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"Nexus_S2_MyC_SCC.pdf"),width = 4,height = 4)
ggplot(s2_myc,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+ theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+  
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")
dev.off()


## nexus tbp k562 and  venters tbp k562
rl <- qc_ind %>% filter(grepl("venters_TBP",name)) %>%
  select(read_length) %>% unique()

tbp_k562 <- scc %>% filter(grepl("K562_TBP",name) |
                                grepl("venters_TBP",name)) %>%
  mutate(cond = ifelse(grepl("chipnexus",name),"ChIP-nexus","ChIP-exo"),
         repl = ifelse(grepl("rep1",tolower(name)),"Rep-1",
           ifelse(grepl("rep2",tolower(name)),"Rep-2","Rep-3")))

pdf(file = file.path(figs_dir,"TBP_K562_SCC.pdf"),width = 6, height = 4)
ggplot(tbp_k562,aes(shift,cross.corr,colour = repl))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = r,name = "Replicate")+
  facet_grid( . ~ cond)+
  geom_vline(xintercept = rl[,mean(read_length)],linetype = 2)+    
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()


## meijsin gr in different cell lines
rl <- qc_ind %>% filter(grepl("meij",name)) %>%
  select(read_length) %>% unique()

mejsing <- scc %>% filter(grepl("meij",name)) %>%
  mutate(cell = ifelse(grepl("K562",name),"K562",
           ifelse(grepl("U2O",name),"U2OS","IMR90")))

pdf(file = file.path(figs_dir,"Meijsing_GR_SCC.pdf"),width = 4, height = 4)
ggplot(mejsing,aes(shift,cross.corr,colour = cell))+
  geom_line(size = .6)+
  theme_bw()+theme(legend.position = "top")+
  scale_color_manual(values = rev(r),name = "Cell line")+
  geom_vline(xintercept = rl[[1]],linetype = 2)+    
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")                    
dev.off()

## landick biosample 1




## landick biosample 2






