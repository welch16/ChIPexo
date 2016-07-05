


rm(list = ls())

library(parallel)
library(data.table)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)

dir1 <- "data/ChIPexo_QC_runs"
dir2 <- "data/ChIPexo_QC_runs_sig70"

files <- c(list.files(dir1,full.names = TRUE),
           list.files(dir2,full.names = TRUE))

## functions

load_file <- function(x){
  load(x)
  ext_stats}

name_and_join <- function(my_list)
{
  if(is.null(names)){
    names(my_list) <- as.character(1:length(my_list))
  }
  my_list <- mcmapply(
    function(x,y)x[,cond := y],my_list,names(my_list),
    SIMPLIFY = FALSE,mc.cores = 10)
  out <- do.call(rbind,my_list)
  return(out)
}


enrichment_plot1 <- function(stats,incSingle = TRUE)
{
  stats <- stats[,.(cond,ave_reads,cover_rate,f,r)]
  r <- viridis::viridis(1e3, option = "D")
  if(!incSingle){
    stats <- stats[f > 0 & r > 0]
  }
  p <- ggplot(stats,aes(ave_reads,cover_rate))+
    stat_binhex(bins = 70)+facet_grid( ~ cond)+
    scale_x_continuous(limits = c(0,4))+
    scale_y_continuous(limits = c(0,1))+
    theme_bw()+theme(legend.position = "top")+
    scale_fill_gradientn(colours = r,
      trans = 'log10',
      labels=trans_format('log10',math_format(10^.x)) )+
    xlab("Average read coverage (ARC)")+
    ylab("Unique read coverage rate (URCR)")
  return(p)
}

data_files <- files[grepl("TBP",files) & !grepl("Rep1",files)]
data_files <- data_files[-6]


stats <- mclapply(data_files,load_file,mc.cores = 20)
stats <- lapply(stats,function(x)x$stats)


plot1 <- function(stats)
{
  stats <- stats[,.(ave_reads,npos,width)]
  r <- viridis::viridis(1e3, option = "D")
  p <- ggplot(stats,aes(ave_reads,npos /width))+
    stat_binhex(bins = 70)+
    theme_bw()+theme(legend.position = "top")+
    scale_fill_gradientn(colours = r,
      trans = 'log10',
      labels=trans_format('log10',math_format(10^.x)) )+
    xlim(0,100)+ylim(0,2)   
  return(p)
}

plots <- mclapply(stats,plot1,mc.cores = 20)


pdf("test.pdf")
u <- mapply(function(x,y)print(x + ggtitle(y)),plots,basename(data_files),
            SIMPLIFY = FALSE)
dev.off()








