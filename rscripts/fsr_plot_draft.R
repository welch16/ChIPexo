
rm(list = ls())

## the idea of this script is to compare the enriched areas across three samples

library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

folders <- paste0("FoxA1-rep",1:3)
dir <- "/p/keles/ChIPexo/volume3/Analysis/Carroll/mouse"
figs_dir <- "figs/reproducibility/Carroll_mouse"


build_files <- function(folder,dir)
{
  files <- list.files(file.path(dir,folder,"data"))
  patterns <- c("depth","reads_by_region","regions","summary_stat")
  files <- files[sapply(patterns,function(x,files)
                        grep(x,files),files)]
  return(files)
}

files <- lapply(folders,build_files,dir)
mc <- detectCores()

load_files <- function(files,folder,dir)
{
  dir1 <- file.path(dir,folder,"data")

  load(file.path(dir1,files[1]))  ## depth
  load(file.path(dir1,files[2]))  ## reads_table
  load(file.path(dir1,files[3]))  ## regions
  load(file.path(dir1,files[4]))  ## summary_stats

  out <- list()
  out[["depth"]] <- depth
  out[["reads"]] <- reads_table
  out[["regions"]] <- regions
  out[["summary"]] <- summary_stats

  return(out)
                
}

## load the data
replicates <- mapply(load_files,
  files,folders,MoreArgs = list(dir),
  SIMPLIFY = FALSE)                       


## clean reads tables and set keys
replicates <- lapply(replicates,function(x){
  x[["reads"]] <- lapply(x[["reads"]],
    function(y){
    y[,match := paste(seqnames,region,sep ="_")]
    y[,region := NULL]
    return(y)})
  x[["reads"]] <- do.call(rbind,x[["reads"]])
  setkey(x[["reads"]],match)
  return(x)})
names(replicates) <- folders
chr <- names(replicates[[1]]$regions)


## tidy regions
replicates <- lapply(replicates,
  function(x,chr){
    out <- mclapply(chr,function(y,regions){
      dt <- data.table(seqnames = y ,
                       start = start(regions[[y]]),
                       end = end(regions[[y]]))
      dt[,match := paste(seqnames,1:nrow(dt),sep = "_")]
      return(dt)},x[["regions"]],mc.cores = mc)
    out <- do.call(rbind,out)
    x[["regions"]] <- out
    return(x)},chr)

replicates <- lapply(replicates,function(x){
  out <- mclapply(x[["summary"]],function(y){
    y[,region := paste(chrID,region,sep = "_")]
    y[,chrID := NULL]
    nms <- names(y)
    setnames(y,nms,c("match",nms[-1]))
    return(y)},mc.cores = mc)
  out <- do.call(rbind,out)
  x[["summary"]] <- out
  setkey(x[["summary"]],match)
  return(x)})

replicates <- lapply(replicates,function(x){
  x[["stats"]] <- merge(x[["regions"]],x[["summary"]],by = "match")
  x[["summary"]] <- NULL
  x[["regions"]] <- NULL
  return(x)})

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
    scale_color_brewer(name = "quantiles",palette = "Set1")+ylim(0,1)+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top")+ylab("fwd strand ratio")+xlab("least amount of fragments in region")
  return(p)
}

stats <- lapply(replicates,function(x)x$stats)

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
  dt[ ,variable := factor(variable, levels = ord)]
  p <- ggplot(dt , aes(values,value, fill = variable))+geom_bar(stat="identity")+
    scale_fill_brewer(name = "Strand composition",palette = "Pastel1")+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top")+
    xlab("least amount of fragments in region")
  if(prop){
    p <- p + ylab("Cummulative proportion")+ylim(0,1)
  }else{
    p <- p + ylab("Cummulative nr. regions")+scale_y_log10(label = trans_format('log10',math_format(10^.x)))
  }

  return(p)
}

values <- seq(5,150,by = 5)
plots <- lapply(stats,label_plot,values = values,prop = TRUE)
pdf(file = "figs/fsr/draft_label_prop_plot.pdf",width = 9,height = 4)
u <- lapply(plots,print)
dev.off()


plots <- lapply(stats,label_plot,values = values)
pdf(file = "figs/fsr/draft_label_counts_plot.pdf",width = 9,height = 4)
u <- lapply(plots,print)
dev.off()


plots <- lapply(stats,forward_strand_ratio_plot,probs = c(.05,.95,.75,.25,.5),values = values)
pdf(file = "figs/fsr/draft_fsr_plot.pdf",width = 9,height= 4)
u <- lapply(plots,print)
dev.off()
