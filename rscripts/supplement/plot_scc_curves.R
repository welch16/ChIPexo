
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

rl_exo = 58
rl_set = 32


landick1 <- scc %>% filter(grepl("land",name)) %>%
  filter(grepl("93",name) | grepl("O2",name)) %>%
  mutate(prot = ifelse(grepl("SE",name),"SE ChIP-seq",
           ifelse(grepl("PE",name),"PE ChIP-seq","ChIP-exo")),
         repl = ifelse(grepl("933",name) | grepl("937",name),"Rep-2","Rep-1"),
         cond = ifelse(grepl("70-O2",name),"Exp - O_2",
           ifelse(grepl("935",name) | grepl("937",name),"Stat + O_2","Exp + O_2")))

pdf(file = file.path(figs_dir,"Landick_old_SCC.pdf"),width = 9, height = 4)
ggplot(landick1,aes(shift,cross.corr,colour = cond,linetype = repl))+
  geom_line(size = .6)+theme_bw()+facet_grid(~  prot)+
  scale_color_manual(values = rev(r),name = "Growth condition")+
  scale_linetype_manual(values = c(1,2),name = "Replicate")+
  theme(legend.position = "top")+
  geom_vline(xintercept =c( rl_exo,rl_set),linetype = 2)+    
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")      
dev.off()


## landick biosample 2

rl_exo = 51
rl_set = 51

landick2 <- scc %>% filter(grepl("land",name)) %>%
  filter(grepl("139",name) | grepl("140",name) | grepl("131",name) | grepl("1320",name)) %>%
  mutate(prot = ifelse(grepl("SE",name),"SE ChIP-seq",
           ifelse(grepl("PE",name),"PE ChIP-seq","ChIP-exo")),
         repl = ifelse(grepl("1311",name) | grepl("1314",name) |
           grepl("1396",name) | grepl("1398",name),"Rep-1","Rep-2"),
         rif = ifelse(grepl("1311",name) | grepl("1317",name) |
           grepl("1396",name) | grepl("1400",name),"0 min.","20 min."))
         
                                                                           
pdf(file = file.path(figs_dir,"Landick_new_SCC.pdf"),width = 9, height = 4)
ggplot(landick2,aes(shift,cross.corr,colour = rif,linetype = repl))+
  geom_line(size = .6)+theme_bw()+facet_grid(~  prot)+
  scale_color_manual(values = rev(r),name = "Rif. time")+
  scale_linetype_manual(values = c(1,2),name = "Replicate")+
  theme(legend.position = "top")+
  geom_vline(xintercept =c( rl_exo,rl_set),linetype = 2)+    
  xlab("Shift")+ylab("Strand Cross-Correlation (SCC)")      
dev.off()
         


