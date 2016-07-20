
rm(list = ls())

library(dplyr)
library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)
library(ggplot2)
library(scales)
library(RColorBrewer)

## parameters for plots
window_length <- 25
sm <- 1
topM <- 5000

## fimo results

fimo_dr <- "/p/keles/ChIPexo/volume4/tbp_analysis/fimo"
files <- list.files(fimo_dr,pattern = "fimo.txt",recursive = TRUE,
                    full.names = TRUE)

fimo <- lapply(files,fread)
names(fimo) <-  gsub("_jaspar","",basename(dirname(files)))

fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("motifID","sequenceID",
                        "motifStart","motifEnd","strand",
                        "score","pval","qval","sequence"))
  return(x)})

chr_from_fimo <- function(x)sapply(strsplit(as.character(x),":"),
                function(y)y[1])

start_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[1])
  return(as.numeric(out))
}

end_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[2])
  return(as.numeric(out))
}

fimo <- lapply(fimo,function(x){
  x[,seqnames := chr_from_fimo(sequenceID)]
  x[,start := start_from_fimo(sequenceID)]
  x[,end := end_from_fimo(sequenceID)]
  return(x)})

## topM

mm <- c(50,100,250,500,1000,2000,4000,8000)

get_scores <- function(fimo,topM)fimo[order(-score)][1:topM]

summarize <- function(topM,fimo){
  scores <- lapply(fimo,get_scores,topM)
  scores <- mapply(function(x,y)x[,cond := y],scores,
                   names(scores),SIMPLIFY = FALSE)
  scores <- do.call(rbind,scores)
  return(scores)
}

scores <- mclapply(mm,summarize,fimo,mc.cores = 20)
scores <- mapply(function(x,y)x[,topM := y],scores,
                 mm,SIMPLIFY = FALSE)
scores <- do.call(rbind,scores)

scores <- scores %>%
  mutate(prot = ifelse(grepl("venters",cond),"ChIP-exo","ChIP-nexus"),
         repl = ifelse(grepl("Rep1",cond),"Rep-1",
           ifelse(grepl("Rep2",cond),"Rep-2","Rep-3")),
         topM = factor(topM))

scores <- data.table(scores)

r <- brewer.pal(n= 8,"Set1")



pdf(file = "figs/for_paper/TBP_scores_comparison.pdf",height = 4,width = 8)
ggplot(scores,aes(cond,score,colour = prot))+
  geom_boxplot(outlier.size = 0.1)+facet_grid( . ~ topM)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5),
        axis.title.x = element_blank(),
        ## strip.background = element_blank(),
        legend.position = "top")+
  ## scale_x_discrete(labels = paste("Rep",c(1:2,1:3),sep = "-"))+
  scale_color_manual(name = "",values = r[c(4,7)])+
  ylab("FIMO score")
dev.off()
