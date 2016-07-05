
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
    xlab("Average Read Coefficient (ARC)")+
    ylab("Unique Read Coefficient (URC)")
  return(p)
}

fsr_plot1 <- function(stats,probs = c(0,.25,.5,.75,1),
                      values = 1:750,incSingle = TRUE)
{  
  if(!incSingle){
    stats <- stats[f > 0 & r > 0]
  }
  stats <- stats[,.(cond,depth,fsr)]
  rows <- mclapply(values,
    function(i)stats[depth > i,quantile(fsr,probs),by = cond],mc.cores = 20)
  rows <- mcmapply(function(rr,i,quant){
    rr <- rr[, depth := i]
    nn <- rr[,length(unique(cond))]
    rr <- rr[,quantile := rep(quant,nn)]
    return(rr)},rows,values,MoreArgs = list(probs),SIMPLIFY = FALSE,
    mc.cores = 20)
  rows <- do.call(rbind,rows)
  setnames(rows,names(rows),c("cond","fsr","depth","quantiles"))
  out <- ggplot(rows,aes(depth,fsr,colour = as.factor(quantiles)))+
    geom_line(size = 1)+
    theme_bw()+theme(legend.position = "top")+facet_grid(cond ~ .)+
    scale_color_brewer(palette = "Dark2",name = "")+
    xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
    scale_x_continuous(limits = c(0,300),breaks = seq(0,300,by = 50))+
    ylim(0,1)
  return(out)
}

label_plot1 <- function(stats,values = 1:750)
{
  stats <- stats[,.(depth,f,r,cond)]
  stats[,label := "both"]
  stats[,label := ifelse(f > 0 & r > 0 , label,
           ifelse( f > 0 ,"fwd","bwd"))]
  setkey(stats,label)
  fwd <- stats["fwd"]
  fwd <- mclapply(values,
    function(i)fwd[depth > i , length(depth),by = cond],
    mc.cores = 20)
  bwd <- stats["bwd"]
  bwd <- mclapply(values,
    function(i)bwd[depth > i , length(depth),by = cond],
    mc.cores = 20)
  both <- stats["both"]
  both <- mclapply(values,
    function(i)both[depth > i , length(depth),by = cond],
    mc.cores = 20)
  dt <- mapply(function(val,fwd,bwd,both){
    dt <- rbind(fwd[,label := "fwd"],
                     bwd[,label := "bwd"],
                     both[,label := "both"])
    sums <- dt[,sum(V1),by = cond]
    dt <- merge(dt,sums,by = "cond")
    dt[,prop := V1.x / V1.y]
    dt[,depth := val]
    dt[,V1.x := NULL]
    dt[,V1.y := NULL]    
    return(dt)},values,fwd,bwd,both,SIMPLIFY = FALSE)
  dt <- do.call(rbind,dt)

  r <- brewer.pal(name = "Set1",3)
  dt[,label := factor(label,levels = c("both","fwd","bwd"))]
  out <- ggplot(arrange(dt,label, desc(prop)),
                aes(depth,prop,fill = label))+
    geom_bar(stat = "identity")+
    theme_bw()+theme(legend.position = "top")+
    facet_grid(cond ~ .)+
    scale_fill_brewer(palette = "Pastel1",name = "Strand composition")+
    xlab("Minimum number of reads")+ylab("Proportion of regions")
  return(out)
}

figsdir <- "figs/supplement"

values <- 1:300
probs <- c(.1,.9,.75,.25,.5)

bw1 <- 3
bh1 <- 4

bw2 <- 6
bh2 <- 1.5


## pugh ctcf
prefix <- "Pugh_CTCF_hela"
data_files <- files[grep("pugh",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("CTCF")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 1 ,height = bh1)
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 1.4)
u <- lapply(plots[3:5],print)
dev.off()


## carroll human
prefix <- "Carroll_ER_MCF7"
data_files <- files[grep("carroll_human",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2","Rep-3")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 3 ,height = bh1)
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 3)        
u <- lapply(plots[3:5],print)
dev.off()

## carroll mouse

prefix <- "Carroll_FoxA1_mouse"
data_files <- files[grep("carroll_mouse",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2","Rep-3")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
          width =  bw1 * 3, height = bh1)
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 3)      
u <- lapply(plots[3:5],print)
dev.off()


## nexus embryo dorsal

prefix <- "ChIPnexus_embryo_dorsal"
data_files <- files[grep("embryo_dorsal",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")), 
      width =  bw1 * 2 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
            width = bw2 ,height = bh2 * 2)
u <- lapply(plots[3:5],print)
dev.off()

## nexus embryo twist

prefix <- "ChIPnexus_embryo_twist"
data_files <- files[grep("embryo_twist",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 2 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 2)    
u <- lapply(plots[3:5],print)
dev.off()


## nexus s2 max

prefix <- "ChIPnexus_S2_max"
data_files <- files[grep("S2_Max",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 2 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 2)    
u <- lapply(plots[3:5],print)
dev.off()


## nexus s2 myc

prefix <- "ChIPnexus_S2_myc"
data_files <- files[grep("S2_MyC",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 2 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 2)  
u <- lapply(plots[3:5],print)
dev.off()


## nexus tbp k562

prefix <- "ChIPnexus_K562_TBP"
data_files <- files[grepl("K562_TBP",files) & !grepl("sample",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 2 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 2)  
u <- lapply(plots[3:5],print)
dev.off()

## venters tbp k562

prefix <- "Venters_K562_TBP"
data_files <- files[grepl("TBP_K562",files) & !grepl("sample",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("Rep-1","Rep-2","Rep-3")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 3 ,height = bh1)    
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 3)    
u <- lapply(plots[3:5],print)
dev.off()

## meijsin gr in different cell lines

prefix <- "Meijsing_GR"
data_files <- files[grepl("meij",files)]
data <- mclapply(data_files,load_file,mc.cores = 3)
data <- lapply(data,function(x)x$stats)
names(data) <- c("IMR90","K562","U2OS")
data <- name_and_join(data)
plots <- list()
plots[[1]] <- enrichment_plot1(data,TRUE)
plots[[2]] <- enrichment_plot1(data,FALSE)
plots[[3]] <- fsr_plot1(data,probs,values,TRUE)
plots[[4]] <- fsr_plot1(data,probs,values,FALSE)
plots[[5]] <- label_plot1(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width =  bw1 * 3 ,height = bh1)         
u <- lapply(plots[1:2],print)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = bw2 ,height = bh2 * 3)  
u <- lapply(plots[3:5],print)
dev.off()

