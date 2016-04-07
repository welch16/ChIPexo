rm(list = ls())
library(viridis)
library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(hexbin)
library(reshape2)

load_all("~/Desktop/Docs/Code/ChIPexoQual")


indir <- "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert= TRUE)]
mc <- detectCores()


files <- files[grep("TBP",files)]
prefix <- "Venters_TBP_human"

expts <- lapply(file.path(indir,files),create_exo_experiment,calc_summary = TRUE,height = 1,parallel = TRUE , mc = detectCores())

names(expts) <- gsub(".sort.bam","",files)
stats <- lapply(expts,summary_stats)


fsr_DT <- function(stats,probs = c(0,.25,.5,.75,1),values = 1:750, mc = 8)
{
  stopifnot(is.numeric(values))
  stopifnot(length(probs) <= 9)
  stopifnot(length(values) > 1)
  to_use <- stats[,.(depth, fsr)]
  rows <- mclapply(values,function(i)to_use[depth > i,quantile(fsr, probs)],mc.cores = mc)
  rows <- mapply(function(x,i)data.table(depth = i,quantiles = probs,fsr = x ),rows,values,SIMPLIFY = FALSE)
  rows <- do.call(rbind,rows)
  rows[,quantiles := factor(quantiles,levels = sort(probs))]
  return(rows)
}


label_DT <- function(stats,values = 1:750,mc = 8,prop = FALSE)
{
  stopifnot(is.numeric(values))
  stopifnot(length(values) > 1)

  stats[,label := ifelse( f > 0 & r > 0 , "both",ifelse(f > 0,"fwd","bwd"))]

  to_use <- stats[,.(depth,label)]
  setkey(to_use,label)
  fwd <- to_use["fwd"]
  fwd <- do.call(c,mclapply(values,function(i)nrow(fwd[depth > i]),mc.cores = mc))
  bwd <- to_use["bwd"]
  bwd <- do.call(c,mclapply(values,function(i)nrow(bwd[depth > i]),mc.cores = mc))
  both <- to_use["both"]
  both <- do.call(c,mclapply(values,function(i)nrow(both[depth > i]),mc.cores = mc))
  dt <- data.table(fwd,bwd,both)  
  if(prop){
    rs <- rowSums(dt)
    dt <- dt / rs
  }
  ord <- c("both","fwd","bwd")
  setcolorder(dt,ord)
  dt <- cbind(values,dt)
  dt <- melt(dt,id.vars = "values")
  dt[ ,variable := factor(variable, levels = ord)]
  return(dt)
}


values <- seq(1,300,by = 1)
fsr <- lapply(stats,fsr_DT,probs = c(.1,.9,.75,.25,.5),values = values,mc = mc)

stats2 <- lapply(stats,function(x)x[f > 0 & r > 0])
fsr2 <- lapply(stats2,fsr_DT,probs = c(.1,.9,.75,.25,.5),values = values,mc = mc)

values <- 1:50
label <- lapply(stats,label_DT,values = values,prop = TRUE)

name_and_join <- function(list)
{
  if(is.null(names)){
    names(list) <- as.character(1:length(list))
  }

  list <- mapply(function(x,y)x[,cond := y],list,names(list),SIMPLIFY = FALSE)
  out <- do.call(rbind,list)

  return(out)
}

stats <- name_and_join(stats)
fsr <- name_and_join(fsr)
label <- name_and_join(label)
fsr2 <- name_and_join(fsr2)
stats2 <- name_and_join(stats2)

figsdir <- "figs/ChIPexoQuals"

r <- viridis::viridis(1e3, option = "D")

pdf(file = file.path(figsdir,paste(prefix,"enrichment.pdf",sep = "_")),width =9 ,height = 4)
ggplot(stats,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+facet_grid( ~ cond)+
  scale_x_continuous(limits = c(0,4))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"enrichment_noSingleStrand.pdf",sep = "_")),width =9 ,height = 4)
ggplot(stats2,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+facet_grid( ~ cond)+
  scale_x_continuous(limits = c(0,4))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"MA_plot.pdf",sep = "_")),width =9 ,height = 4)
ggplot(stats2,aes(M,A))+stat_binhex(bins = 70)+facet_grid( ~ cond)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+
  xlab("M")+ylab("A")+ylim(-10,10)+xlim(-15,5)
dev.off()


pdf(file = file.path(figsdir,paste(prefix,"fsr_surv.pdf",sep ="_")))
ggplot(fsr,aes(depth,fsr,colour = as.factor(quantiles)))+geom_line(size = 1)+
  theme_bw()+theme(legend.position = "top")+facet_grid(cond ~ .)+
  scale_color_brewer(palette = "Dark2",name = "")+
  xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
  xlim(0,300)+ylim(0,1)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"fsr_surv_noSingleStrand.pdf",sep ="_")))
ggplot(fsr2,aes(depth,fsr,colour = as.factor(quantiles)))+geom_line(size = 1)+
  theme_bw()+theme(legend.position = "top")+facet_grid(cond ~ .)+
  scale_color_brewer(palette = "Dark2",name = "")+
  xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
  xlim(0,300)+ylim(0,1)
dev.off()

pdf(file = file.path(figsdir,paste(prefix,"label_surv.pdf",sep ="_")))
ggplot(label,aes(values,value,fill = variable))+geom_bar(stat = "identity")+
  theme_bw()+theme(legend.position = "top")+facet_grid(cond ~ .)+
  scale_fill_brewer(palette = "Pastel1",name = "Strand composition")+
  xlab("Minimum number of reads")+ylab("Proportion of islands")
dev.off()


