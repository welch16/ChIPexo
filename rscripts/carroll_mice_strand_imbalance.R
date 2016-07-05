
rm(list = ls())

library(reshape2)
library(ChIPUtils)
library(ggplot2)
library(data.table)
library(GenomicAlignments)
library(parallel)
library(viridis)

dr <- "data/ChIPexo_QC_runs"
files <- list.files(dr,full.names = TRUE)
files <- files[grepl("carr",files) & grepl("mouse",files)]

load_stat <- function(x){
  load(x)
  ext_stats$stats}

stats <- mclapply(files,load_stat,mc.cores = 3)



forward_strand_ratio_plot <- function(stat,probs = c(0,.25,.5,.75,1),values = 1:750, mc = 8)
{
  stopifnot(is.numeric(values))
  stopifnot(length(probs) <= 9)
  stopifnot(length(values) > 1)
  to_use <- stat[,.(depth, fsr)]
  rows <- mclapply(values,function(i)to_use[depth > i,quantile(fsr, probs)],mc.cores = mc)
  rows <- mapply(function(x,i)data.table(depth = i,quantiles = probs,fsr = x ),rows,values,SIMPLIFY = FALSE)
  rows <- do.call(rbind,rows)
  rows[,quantiles := factor(quantiles,levels = sort(probs))]
  return(rows)
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
  dt[ ,variable := factor(variable, levels = ord)]
  return(dt)

}

names(stats) <- c("Rep-1","Rep-2","Rep-3")

fsr_dt <- mclapply(stats,forward_strand_ratio_plot,probs = c(.1,.25,.5,.75,.9),values = 1:300,mc.cores =3)

label_dt <- lapply(stats,label_plot,prop = TRUE,values = 1:50)

fsr_dt <- mapply(function(x,y)x[,sample := y],fsr_dt,names(fsr_dt),SIMPLIFY = FALSE)

label_dt <- mapply(function(x,y)x[,sample := y],label_dt,names(label_dt),SIMPLIFY = FALSE)

fsr_dt <- do.call(rbind,fsr_dt)
label_dt <- do.call(rbind,label_dt)


figs_dir <- "figs/Carroll_mice_for_paper/"


library(grid)
library(gridExtra)


## fsr_dt[, sample := plyr::mapvalues(sample,
##     from = c("ERR336935.bam","ERR336942.bam","ERR336956.bam"),
##     to = c("rep-3","rep-1","rep-2"))]

label_dt[, sample := plyr::mapvalues(sample,
    from = c("ERR336935.bam","ERR336942.bam","ERR336956.bam"),
    to = c("rep-3","rep-1","rep-2"))]



pdf(file = file.path(figs_dir,"Strand_imbalance.pdf"),width = 9 , height = 5)
p1 <- ggplot(fsr_dt,aes(depth , fsr , colour = quantiles))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1")+ylim(0,1)+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  ylab("Fwd strand ratio (FSR)")+xlab("Min number of reads")+
  facet_grid( sample ~ .)+ggtitle("A")
p2 <- ggplot(label_dt, aes(values,value, fill = variable))+geom_bar(stat="identity")+
  scale_fill_brewer(name = "",palette = "Pastel1")+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+facet_grid(sample ~ .) +
  xlab("Min number of reads")+ggtitle("B")+ylab("Proportion of islands")
print(p1)
print(p2)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1,1:2))
print(p2, vp = vplayout(1,3))
dev.off()


du <- data.table(bam = basename(files))
du[,peaks := c("do2435","do2160","do2161")]

peak_dir <- "/p/keles/ChIPexo/volume4/carroll_data/peaks"

peak_files <- list.files(peak_dir)

peak_files <- peak_files[sapply(du[,(peaks)],grep,peak_files)]

peaks <- lapply(file.path(peak_dir,peak_files) ,read.table)

peaks <- lapply(peaks,data.table)
names(peaks) <- du[,(peaks)]

peak_stat_ov <- function(peaks,stats,nm,minDep  =1)
{
  peaks <- copy(peaks[,1:3,with = FALSE])
  setnames(peaks,names(peaks),c("seqnames","start","end"))
  peaks <- dt2gr(peaks)
  stats <- copy(stats[depth > minDep, .(seqnames,start,end,fsr)])
  stats[,strand := "*"]
  gr <- dt2gr(stats)
  mcols(gr)$ov <- ifelse(countOverlaps(gr,peaks) > 0 ,1 ,0)
  dt <- data.table(as.data.frame(mcols(gr)))
  dt[, ov := ifelse(ov > 0 , "yes","no")]
  dt[,sample := nm]
  test <- ks.test(dt[ov == "yes",(fsr)],dt[ov == "no",(fsr)])
  out <- list(dt,test)
  return(out)
}

low <- mcmapply(peak_stat_ov,peaks,stats,gsub(".bam","",du[,(bam)]),MoreArgs = list(minDep = 1),mc.cores = 3,SIMPLIFY = FALSE)
med <- mcmapply(peak_stat_ov,peaks,stats,gsub(".bam","",du[,(bam)]),MoreArgs = list(minDep = 30),mc.cores = 3,SIMPLIFY = FALSE)
high <- mcmapply(peak_stat_ov,peaks,stats,gsub(".bam","",du[,(bam)]),MoreArgs = list(minDep = 100),mc.cores = 3,SIMPLIFY = FALSE)


