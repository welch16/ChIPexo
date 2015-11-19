
rm(list = ls())

library(reshape2)
library(ChIPUtils)
library(ggplot2)
library(data.table)
library(GenomicAlignments)
library(parallel)
library(viridis)

dr <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(dr)

files <- files[grep("bam",files)]
files <- files[grep("bai",files,invert = TRUE)]

files <- file.path(dr,files)

reads <- mclapply(files,create_reads,mc.cores = 3)
names(reads) <- basename(files)


create_regions <- function(reads,lower,rangesOnly=TRUE)
{
  stopifnot(class(reads) == "GRanges")
  stopifnot(lower > 0)
  sreads <- split(reads,seqnames(reads))
  sreads <- as.list(sreads)
  cover <- mcmapply(function(x,y)coverage(x)[[y]],sreads,names(sreads),mc.cores = detectCores())
  islands <- lapply(cover,slice,lower = lower, rangesOnly = rangesOnly)
  islands <- mapply(function(x,y)GRanges(seqnames = y,ranges = x),islands,names(islands),SIMPLIFY =FALSE)
  islands <- do.call(rbind,lapply(islands,gr2dt))
  return(islands)
}

gr <- lapply(reads,function(x){
  byChr <- mapply(rbind,readsF(x),readsR(x),SIMPLIFY =FALSE)
  out <- do.call(rbind,byChr)
  return(dt2gr(out))
})

regs <- lapply(gr,create_regions,lower  = 1)

build_stats <- function(region,reads)
{
  ## fix formats and stuff
  region <- dt2gr(region)
  
  ov <- findOverlaps(region,reads)
  reads <- gr2dt(reads)
  w <- width(region)    
  region <- gr2dt(region)
  region[ , width := w]
  region[, match := paste0(seqnames,":",start,"_",end)]
  reads[  subjectHits(ov), match := region[queryHits(ov), (match)] ]  
  reads[,strand := ifelse(strand == "+", "F","R")]
  reads <- reads[!is.na(match)]

  ## get base statistics
  f <- reads[,sum(strand == "F"),by = match]
  setnames(f,names(f),c("match","f"))
  setkey(f,match)
  r <- reads[,sum(strand == "R"),by = match]
  setkey(r,match)
  setnames(r,names(r),c("match","r"))
  f_uniq <- reads[strand == "F",length(unique(start)),by = match]
  setnames(f_uniq,names(f_uniq),c("match","f_pos"))
  setkey(f_uniq,match)
  r_uniq <- reads[strand == "R",length(unique(end)),by = match]
  setnames(r_uniq,names(r_uniq),c("match","r_pos"))
  setkey(r_uniq,match)

  ## merge statistics
  stats <- merge(region,f,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,r,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,f_uniq,by = "match",allow.cartesian = TRUE,all = TRUE)
  stats <- merge(stats,r_uniq,by = "match",allow.cartesian = TRUE, all = TRUE)
  stats[is.na(f_pos), f_pos := 0]
  stats[is.na(r_pos), r_pos := 0]

  ## calculate composite stats
  stats[ , depth := f + r]
  stats[ , npos := f_pos + r_pos]
  stats[ , ave_reads := depth / width]
  stats[ , cover_rate := npos / depth]
  stats[ , fsr := f / (f + r)]
  stats[ , label := ""]
  stats[ f > 0 & r > 0 , label := "both"]
  stats[ f > 0 & r == 0, label := "fwd"]
  stats[ r > 0 & f == 0, label := "bwd"]

  stats[ , M := as.numeric(NA)]
  stats[ , A := as.numeric(NA)]

  stats[f > 0 & r > 0, M := log2(f * r) - 2 * log2(width)]
  stats[f > 0 & r > 0, A := log2( f/ r)]

  stats[ , strand := NULL]

  return(stats)
}

stats <- mcmapply(build_stats,regs,gr,mc.cores = 6,SIMPLIFY = FALSE)


## common <- subsetByOverlaps(dt2gr(regs[[1]]),subsetByOverlaps(dt2gr(regs[[2]]),dt2gr(regs[[3]])))

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

fsr_dt <- mclapply(stats,forward_strand_ratio_plot,probs = c(.1,.25,.5,.75,.9),values = 1:300,mc.cores =3)
label_dt <- lapply(stats,label_plot,prop = TRUE,values = 1:50)

fsr_dt <- mapply(function(x,y)x[,sample := y],fsr_dt,names(fsr_dt),SIMPLIFY = FALSE)
label_dt <- mapply(function(x,y)x[,sample := y],label_dt,names(label_dt),SIMPLIFY = FALSE)

fsr_dt <- do.call(rbind,fsr_dt)
label_dt <- do.call(rbind,label_dt)


figs_dir <- "figs/Carroll_mice_for_paper/"


library(grid)
library(gridExtra)



pdf(file = file.path(figs_dir,"Strand_imbalance.pdf"),width = 9 , height = 5)
p1 <- ggplot(fsr_dt,aes(depth , fsr , colour = quantiles))+geom_line(size = 1)+
  scale_color_brewer(name = "",palette = "Set1")+ylim(0,1)+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  ylab("fwd strand ratio")+xlab("least amount of fragments in region")+
  facet_grid( sample ~ .)+ggtitle("A")
p2 <- ggplot(label_dt, aes(values,value, fill = variable))+geom_bar(stat="identity")+
  scale_fill_brewer(name = "",palette = "Pastel1")+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+facet_grid(sample ~ .) +
  xlab("least amount of fragments in region")+ggtitle("B")
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

peak_dir <- "/p/keles/ChIPexo/volume3/CarrollData/peaks"
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


