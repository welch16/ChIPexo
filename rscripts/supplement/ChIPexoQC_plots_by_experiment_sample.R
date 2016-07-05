
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

enrichment_plot2 <- function(stats,incSingle = TRUE,xl = 4)
{
  stats <- stats[,.(cond,samp,ave_reads,cover_rate,f,r)]
  r <- viridis::viridis(1e3, option = "D")
  if(!incSingle){
    stats <- stats[f > 0 & r > 0]
  }
  p <- ggplot(stats,aes(ave_reads,cover_rate))+
    stat_binhex(bins = 70)+facet_grid(samp ~ cond)+
    scale_x_continuous(limits = c(0,xl))+
    scale_y_continuous(limits = c(0,1))+
    theme_bw()+theme(legend.position = "top")+
    scale_fill_gradientn(colours = r,
      trans = 'log10',
      labels=trans_format('log10',math_format(10^.x)) )+
    xlab("Average Read Coefficient (ARC)")+
    ylab("Unique Read Coefficient (URC)")
  return(p)
}

fsr_plot2 <- function(stats,probs = c(0,.25,.5,.75,1),
                      values = 1:750,incSingle = TRUE)
{
  if(!incSingle){
    stats <- stats[f > 0 & r > 0]
  }
  stats <- stats[,.(cond,samp,depth,fsr)]
  rows <- mclapply(values,
    function(i)stats[depth > i,quantile(fsr,probs),by = .(samp,cond)],
    mc.cores = 20)
  rows <- mcmapply(function(rr,i,quant){
    rr <- rr[, depth := i]
    nn <- rr[,length(unique(paste0(cond,samp)))]
    rr <- rr[,quantile := rep(quant,nn)]
    return(rr)},rows,values,MoreArgs = list(probs),SIMPLIFY = FALSE,
    mc.cores = 20)
  rows <- do.call(rbind,rows)
  setnames(rows,names(rows),c("samp","cond","fsr","depth","quantiles"))
  out <- ggplot(rows,aes(depth,fsr,colour = as.factor(quantiles)))+
    geom_line(size = 1)+
    theme_bw()+theme(legend.position = "top")+facet_grid(samp ~ cond)+
    scale_color_brewer(palette = "Dark2",name = "")+
    xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
    xlim(0,max(values))+ylim(0,1)
  return(out)
}

label_plot2 <- function(stats,values = 1:75)
{
  stats <- stats[,.(depth,f,r,cond,samp)]
  stats[,label := "both"]
  stats[,label := ifelse(f > 0 & r > 0 , label,
           ifelse( f > 0 ,"fwd","bwd"))]
  setkey(stats,label)
  fwd <- stats["fwd"]
  fwd <- mclapply(values,
    function(i)fwd[depth > i , length(depth),by = .(samp,cond)],
    mc.cores = 20)
  bwd <- stats["bwd"]
  bwd <- mclapply(values,
    function(i)bwd[depth > i , length(depth),by = .(samp,cond)],
    mc.cores = 20)
  both <- stats["both"]
  both <- mclapply(values,
    function(i)both[depth > i , length(depth),by = .(samp,cond)],
    mc.cores = 20)
  dt <- mapply(function(val,fwd,bwd,both){
    dt <- rbind(fwd[,label := "fwd"],
                     bwd[,label := "bwd"],
                     both[,label := "both"])
    sums <- dt[,sum(V1),by = .(samp,cond)]
    sums[,aux := paste(samp,cond)]
    sums[,samp := NULL]
    sums[,cond := NULL]
    dt[,aux := paste(samp,cond)]
    dt <- merge(dt,sums,by = "aux")
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
    facet_grid(samp ~ cond)+xlim(0,max(values))+
    scale_fill_brewer(palette = "Pastel1",name = "Strand composition")+
    xlab("Minimum number of reads")+ylab("Proportion of regions")
  return(out)
}

figsdir <- "figs/supplement"

## TBP sample
values <- 1:300
probs <- c(.1,.9,.75,.25,.5)

prefix <- "TBP_sample"
data_files <- files[grepl("TBP",files) & grepl("samp",files)]
data_files <- data_files[-5]
#data_files <- data_files[!grepl("Rep1",data_files)]

data <- mclapply(data_files,load_file,mc.cores = 20)
data <- lapply(data,function(x)x$stats)

data <- mcmapply(function(x,y)x[,samp := y],
  data,rep(paste0(seq(20,50,by = 10),"M"),4),SIMPLIFY = FALSE,mc.cores = 20)

data <- mcmapply(function(x,y)x[,cond := y],
  data,
  paste(c(rep("ChIP-nexus",4),rep("ChIP-exo",12)),
        c(rep("Rep-2",4),rep("Rep-1",4),rep("Rep-2",4),rep("Rep-3",4)),
        sep = " "),
  SIMPLIFY = FALSE,mc.cores = 20)

data <- do.call(rbind,data)


plots <- list()
plots[[1]] <- enrichment_plot2(data,TRUE)
plots[[2]] <- enrichment_plot2(data,FALSE)
plots[[3]] <- fsr_plot2(data,probs,values,TRUE)
plots[[4]] <- fsr_plot2(data,probs,values,FALSE)
plots[[5]] <- label_plot2(data,1:50)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width = 12 ,height = 8)
u <- lapply(plots[1:2],print)
dev.off()
 
pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = 7,height = 6)
u <- lapply(plots[3:5],print)
dev.off()


## sig70 sample

values <- 1:120


prefix <- "Sig70_sample_bios1"
data_files <- files[grepl("Sig70",files) & grepl("samp",files) &
                    grepl("edsn9",files)]

data <- mclapply(data_files,load_file,mc.cores = 20)
data <- lapply(data,function(x)x$stats)

samp <- sapply(strsplit(gsub(".RData","",basename(data_files)),"subsample"),
  function(x)x[2])               
samp <- paste0(as.numeric(samp) / 1e3 , "K")

data <- mcmapply(function(x,y)x[,samp := y],
  data,samp,SIMPLIFY = FALSE,mc.cores = 20)


edsn <- sapply(strsplit(basename(data_files),"_"),function(x)x[1])
cond <- plyr::mapvalues(edsn,
       from = c("edsn931","edsn933","edsn935","edsn937"),
       to = c("Aerobic Rep-1","Aerobic Rep-2",
         "Anaerobic Rep-1","Anaerobic Rep-2"))
         
data <- mcmapply(function(x,y)x[,cond := y],
  data,cond,
  SIMPLIFY = FALSE,mc.cores = 20)

data <- do.call(rbind,data)

data[,samp := factor(samp,levels = paste0(seq(10,100,by = 10),"K"))]


plots <- list()
plots[[1]] <- enrichment_plot2(data,TRUE,xl = 2)
plots[[2]] <- enrichment_plot2(data,FALSE,xl = 2)
plots[[3]] <- fsr_plot2(data,probs,values,TRUE)
plots[[4]] <- fsr_plot2(data,probs,values,FALSE)
plots[[5]] <- label_plot2(data,1:30)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width = 22 ,height = 10)
u <- lapply(plots[1:2],print)
dev.off()
 
pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = 10,height = 9)
u <- lapply(plots[3:5],print)
dev.off()

## sig70 sample

values <- 1:120

prefix <- "Sig70_sample_bios2"
data_files <- files[grepl("Sig70",files) & grepl("samp",files) &
                    !grepl("edsn9",files)]

data <- mclapply(data_files,load_file,mc.cores = 20)
data <- lapply(data,function(x)x$stats)

samp <- sapply(strsplit(gsub(".RData","",basename(data_files)),"subsample"),
  function(x)x[2])               
samp <- paste0(as.numeric(samp) / 1e3 , "K")

data <- mcmapply(function(x,y)x[,samp := y],
  data,samp,SIMPLIFY = FALSE,mc.cores = 20)

edsn <- sapply(strsplit(basename(data_files),"_"),function(x)x[1])
cond <- plyr::mapvalues(edsn,
       from = c("edsn1311","edsn1314","edsn1317","edsn1320"),
       to = c("No Rif. Rep-1","Rif. 20 min Rep-1",
              "No Rif. Rep-2","Rif. 20 min Rep-2"))
         
data <- mcmapply(function(x,y)x[,cond := y],
  data,cond,
  SIMPLIFY = FALSE,mc.cores = 20)

data <- do.call(rbind,data)

data[,samp := factor(samp,levels = paste0(seq(10,100,by = 10),"K"))]
data[,cond := factor(cond,levels = c("No Rif. Rep-1","Rif. 20 min Rep-1",
              "No Rif. Rep-2","Rif. 20 min Rep-2"))]

plots <- list()
plots[[1]] <- enrichment_plot2(data,TRUE,xl = 2)
plots[[2]] <- enrichment_plot2(data,FALSE,xl = 2)
plots[[3]] <- fsr_plot2(data,probs,values,TRUE)
plots[[4]] <- fsr_plot2(data,probs,values,FALSE)
plots[[5]] <- label_plot2(data,1:30)

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),
      width = 22 ,height = 10)
u <- lapply(plots[1:2],print)
dev.off()
 
pdf(file = file.path(figsdir,paste(prefix,"strand_imbalance.pdf",sep = "_")),
        width = 10,height = 9)
u <- lapply(plots[3:5],print)
dev.off()

