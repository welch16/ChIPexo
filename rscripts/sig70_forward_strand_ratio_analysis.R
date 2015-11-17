
rm(list = ls())

library(ggplot2)
library(data.table)
library(ChIPUtils)
library(reshape2)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)

## load data from pipeline

dr <- "../Analysis/Landick/stat-vs-exp/Sig70-exp-rep2/data"

files <- list.files(dr)

load( file.path(dr, files[grep("summ",files)]))

stats <- mcmapply(function(x,y)x[,chr := y],summary_stats,names(summary_stats),SIMPLIFY = FALSE,mc.cores = 20)
stats <- do.call(rbind,stats)

forward_strand_ratio_plot <- function(stat,probs = c(0,.25,.5,.75,1),values = 1:750, mc = 8)
{
  stopifnot(is.numeric(values))
  stopifnot(length(probs) <= 9)
  stopifnot(length(values) > 1)
  to_use <- stat[,.(depth, prob)]
  rows <- mclapply(values,function(i)to_use[depth > i,quantile(prob, probs)],mc.cores = mc)
  rows <- mapply(function(x,i)data.table(depth = i,quantiles = probs,fsr = x ),rows,values,SIMPLIFY = FALSE)
  rows <- do.call(rbind,rows)
  rows[,quantiles := factor(quantiles,levels = sort(probs))]
  p <- ggplot(rows,aes(depth , fsr , colour = quantiles))+geom_line(size = 1)+
    scale_color_brewer(name = "quantiles",palette = "Dark2")+ylim(0,1)+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top",plot.title = element_text(hjust =0))+
    ylab("fwd strand ratio")+xlab("least amount of fragments in region")
  return(p)
}


label_plot <- function(stat,values = 1:750,mc = 8,prop = FALSE)
{
  stopifnot(is.numeric(values))
  stopifnot(length(values) > 1)

  to_use <- stat[,.(depth,label)]
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
  if(!prop){
    dt[value == 0 , value := 1]
  }
  dt[ ,variable := factor(variable, levels = ord)]
  p <- ggplot(dt , aes(values,value, fill = variable))+geom_bar(stat="identity")+
    scale_fill_brewer(name = "Strand composition",palette = "Pastel1")+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top",plot.title = element_text(hjust = 0))+
    xlab("least amount of fragments in region")
  if(prop){
    p <- p + ylab("Cummulative proportion")+ylim(0,1)
  }else{
    p <- p + ylab("Cummulative nr. regions")+scale_y_log10(label = trans_format('log10',math_format(10^.x)))
  }

  return(p)
}


values <- 1:50
quants <- c(.05,.25,.5,.75,.95)

p1 <- forward_strand_ratio_plot(stats,probs =quants,value = values,mc = 12)
p2 <- label_plot(stats,values = values, mc = 12,prop = TRUE)
p3 <- label_plot(stats,values = values, mc = 12,prop = FALSE)

grid.arrange(p1,p3,nrow = 2)
dev.off()



