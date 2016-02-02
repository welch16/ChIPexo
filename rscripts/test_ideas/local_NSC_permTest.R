
rm(list = ls())

library(reshape2)
library(ChIPUtils)
library(ggplot2)
library(data.table)
library(GenomicAlignments)
library(parallel)
library(viridis)
library(grid)
library(gridExtra)

mc <- 24

dr <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(dr)

files <- files[grep("bam",files)]
files <- files[grep("bai",files,invert = TRUE)]

files <- file.path(dr,files)

rep1 <- files[2]
rep2 <- files[3]
rep3 <- files[1]

rep1_reads <- create_reads(rep2)

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

gr <- do.call(rbind,mapply(rbind,readsF(rep1_reads),readsR(rep1_reads),
  SIMPLIFY = FALSE))                           

regions <- create_regions(dt2gr(gr),lower = 1)

regions <- dt2gr(regions)




build_stats <- function(region,reads,mc)
{
  ## fix formats and stuff
  stat_by_chr <- function(reg,read)
  {
     ov <- findOverlaps(reg,read)  
     read <- gr2dt(read)
     w <- width(reg)    
     region <- gr2dt(reg)
     region[ , width := w]
     region[, match := paste0(seqnames,":",start,"_",end)]
     read[  subjectHits(ov), match := region[queryHits(ov), (match)] ]  
     read[,strand := ifelse(strand == "+", "F","R")]
     read <- read[!is.na(match)]

     ## get base statistics
     f <- read[,sum(strand == "F"),by = match]
     setnames(f,names(f),c("match","f"))
     setkey(f,match)
     r <- read[,sum(strand == "R"),by = match]
     setkey(r,match)
     setnames(r,names(r),c("match","r"))
     f_uniq <- read[strand == "F",length(unique(start)),by = match]
     setnames(f_uniq,names(f_uniq),c("match","f_pos"))
     setkey(f_uniq,match)
     r_uniq <- read[strand == "R",length(unique(end)),by = match]
     setnames(r_uniq,names(r_uniq),c("match","r_pos"))
     setkey(r_uniq,match)

     ## merge statistics
     stats <- merge(region,f,by = "match",allow.cartesian = TRUE)
     stats <- merge(stats,r,by = "match",allow.cartesian = TRUE)
     stats <- merge(stats,f_uniq,by = "match",
                    allow.cartesian = TRUE,all = TRUE)
     stats <- merge(stats,r_uniq,by = "match",
                    allow.cartesian = TRUE, all = TRUE)
     stats[is.na(f_pos), f_pos := 0]
     stats[is.na(r_pos), r_pos := 0]

     ## calculate composite stats
     stats[ , depth := f + r]
     stats[ , npos := f_pos + r_pos]
     stats[ , ave_reads := depth / width]
     stats[ , cover_rate := npos / depth]
     stats[ , fsr := f / (f + r)]

     stats[ , M := as.numeric(NA)]
     stats[ , A := as.numeric(NA)]

     stats[f > 0 & r > 0, M := log2(f * r) - 2 * log2(width)]
     stats[f > 0 & r > 0, A := log2( f/ r)]

     stats[ , strand := NULL]

     return(stats)

  }

  reg <- split(region,as.character(seqnames(region)))
  read <- split(reads,as.character(seqnames(reads)))

  stat_chr <- mcmapply(stat_by_chr,reg,read,SIMPLIFY = FALSE,mc.cores = mc)

  stats <- do.call(rbind,stat_chr)

  return(stats)
}

stats <- build_stats(regions,dt2gr(gr),mc = mc)

## only keep regions where local-SCC can be calculated
stats <- stats[ f > 0 & r > 0]
stats <- stats[seqnames != "chrM"]

get_candidate_region <- function(mat,stats)dt2gr(stats[mat,2:4,with = FALSE])

cover_DT <- function(reg,reads,n,perm = FALSE)
{

  gr <- do.call(rbind,mapply(rbind,readsF(reads),readsR(reads),SIMPLIFY =FALSE))
  my_gr <- subsetByOverlaps(dt2gr(gr),reg)
  my_gr <- resize(my_gr,1)

  if(perm)strand(my_gr) <- sample(as.character(strand(my_gr)))

  my_gr <- split(my_gr,as.character(strand(my_gr)))
  fwd <- ranges(my_gr[["+"]])
  bwd <- ranges(my_gr[["-"]])

  fwd <- coverage(fwd)
  bwd <- coverage(bwd)

  fwdDT <- data.table(coord = start(reg):end(reg),tags = 0,strand = "F")
  bwdDT <- data.table(coord = start(reg):end(reg),tags = 0,strand = "R")

  fill_cover <- function(DT,cover)
  {
    z1 <- cumsum(runLength(cover))
    z2 <- c(0,runValue(cover))
    DT[,tags := stepfun(z1,z2)(coord)]
    return(DT)
  }

  fwdDT <- fill_cover(fwdDT,fwd)
  bwdDT <- fill_cover(bwdDT,bwd)

  bwdDT[,tags := -tags]

  DT <- rbind(fwdDT,bwdDT)
  DT[,tags := tags * 1e9 / n]

  return(DT)

}

setkey(stats,match)

candidates <- stats[npos >= 200]
mmm <- candidates[,(match)]

split_vec <- function(d,n)split(d, ceiling(seq_along(d)/n))

regions <- mclapply(mmm,get_candidate_region,stats,mc.cores = mc)

s <- 1:34
scc <- mclapply(regions,function(x)local_strand_cross_corr(rep1_reads,x,shift = s),mc.cores = mc)

scc <- mcmapply(function(x,y)x[,match := y],scc,candidates[,(match)],SIMPLIFY = FALSE)
scc <- do.call(rbind,scc)

covers <- mclapply(regions,cover_DT,rep1_reads,n = nreads(rep1_reads),mc.cores = mc)

covers <- mcmapply(function(x,y)x[,match := y],covers,candidates[,(match)],SIMPLIFY = FALSE,mc.cores = mc)
covers <- do.call(rbind,covers)

setkey(covers,match)
setkey(scc,match)

mmm_split <- split_vec(mmm,20)

figs_dir <- "figs/local_NSC"


pdf(file = file.path(figs_dir,"local_NSC_rep3.pdf"),width = 12,height = 12)
for(mm in mmm_split){
  p <- ggplot(scc[match %in% mm],aes(shift,cross.corr))+
    geom_line()+theme_bw()+theme(legend.position = "none")+geom_vline(xintercept = 34)+
      geom_vline(xintercept = 12,linetype = 2)+facet_wrap( ~ match,nrow = 5,ncol = 4)+ylim(-.2,1)+
      geom_abline(slope =0,intercept =0,linetype =2)
  print(p)
}
dev.off()

pdf(file = file.path(figs_dir,"coverage_rep3.pdf"),width = 14,height = 14)
for(mm in mmm_split){
  p <- ggplot(covers[match %in% mm],aes(coord,tags,colour = strand))+
    geom_step()+theme_bw()+theme(legend.position = "none")+
    scale_color_brewer(palette = "Set1")+facet_wrap( ~ match,nrow = 5,ncol = 4,scales = "free")+
    geom_abline(slope =0,intercept =0,linetype =2)
  print(p)
}
dev.off()


## for(mm in 1:length(mmm_split)){
##   message(mm)  
##   print(candidates[mmm_split[[mm]]])
##   message("-----")
## }

## pdf(file.path(figs_dir,"enrichment.pdf"))
## for(mm in mmm_split){
##   p <- ggplot(candidates[mm],aes(ave_reads,cover_rate,label = match))+geom_point()+
##     geom_text(size = 3.5,hjust = 0,check_overlap = TRUE)+theme_bw()+xlim(0,6)+ylim(0,1)
##   print(p)
## }
## dev.off()



## do this for all k in 2:5
## and calculate aic , bic and dic

## or try k = 2 and check the numbers...

library(mixtools)

models <- list()
for(j in 2:3){
models[[j-1]] <- lapply(1:34,function(k)normalmixEM(scc[shift == k,(cross.corr)], k = j))
}



max_shift <- function(shift,cross.corr)shift[which.max(cross.corr)]

candidates <- merge(candidates,scc[,max_shift(shift,cross.corr),by = match],by = "match")

nms <- names(candidates)
setnames(candidates,nms,c(nms[-length(nms)],"max_shift"))


m <- 10:15

A = normalmixEM(scc[candidates[npos >= 200,(match)]][shift == 11,(cross.corr)],k = 5)
str(A)



pdf(file = file.path(figs_dir,"eval_at_point_rep3.pdf"))
for(k in 1:34){
  p <- ggplot(scc[shift == k],aes(cross.corr))+geom_histogram(bins = 100,aes(y = ..density..),width = .5)+
    stat_density(geom = "line",colour = "red")+ggtitle(paste("shift =",k))
  print(p)
}
dev.off()


pdf(file.path(figs_dir,"eval_all_boxplot_rep3.pdf"))
p <- ggplot(scc,aes(factor(shift,levels = 1:34),cross.corr))+
      stat_boxplot(outlier.shape = NA)+
      geom_abline(slope = 0 , intercept = 0 ,linetype = 2,colour = "red")+
      xlab("shift")+ylim(-.1,.4)
print(p+ggtitle("npos >= 34"))
print( p %+% scc[candidates[npos >= 50,(match)]]+ggtitle("npos >= 50"))
print( p %+% scc[candidates[npos >= 100,(match)]]+ggtitle("npos >= 100"))
print( p %+% scc[candidates[npos >= 200,(match)]]+ggtitle("npos >= 200"))
dev.off()



## pdf(file = file.path(figs_dir,"eval_at_point_npos_gr200.pdf"))
## for(k in 1:34){
##   p <- ggplot(scc[candidates[npos >= 200,(match)]][shift == k],aes(cross.corr))+geom_histogram(bins = 100,aes(y = ..density..),width = .5)+
##     stat_density(geom = "line",colour = "red")+ggtitle(paste("shift =",k))
##   print(p)
## }
## dev.off()


## pdf(file = file.path(figs_dir,"eval_at_point_npos_gr100.pdf"))
## for(k in 1:34){
##   p <- ggplot(scc[candidates[npos >= 100,(match)]][shift == k],aes(cross.corr))+geom_histogram(bins = 100,aes(y = ..density..),width = .5)+
##     stat_density(geom = "line",colour = "red")+ggtitle(paste("shift =",k))
##   print(p)
## }
## dev.off()

BV <- function(cross.corr)sum(abs(diff(cross.corr)))
maxJump <- function(cross.corr)max(abs(diff(cross.corr)))

for(k in 1:5){
  x11()
  print(
ggplot(scc[sample(match,2e3),maxJump(cross.corr),by =match],
       aes(V1))+geom_histogram(bins = 100,aes(y = ..density..),width = .8)+
       stat_density(geom = "line",colour = "red")
)}


ggplot(scc[,maxJump(cross.corr),by =match],
       aes(V1))+geom_histogram(bins = 100,aes(y = ..density..),width = .8)+
       stat_density(geom = "line",colour = "red")




qc <- scc[,.(min(cross.corr),
             max(cross.corr)),by = match]

setnames(qc,names(qc),c("match","min","max","medpos"))

candidates <- merge(candidates,qc,by = "match")
candidates[, range := max - min]



ggplot(candidates , aes(max))+geom_histogram(binwidth = .02,aes(y = ..density..))+
  stat_density(geom = "line",colour = "red")

pdf(file.path(figs_dir,"range_below_readLength.pdf"))


ggplot(candidates[npos >= 200],aes(range))+geom_histogram(binwidth = .02,aes(y = ..density..))+
  stat_density(geom = "line",colour = "red")+xlab("range of local-scc, shift <= 34")


dev.off()





## cc <- local_strand_cross_corr(rep1_reads,reg,shift = 1:150,perm = FALSE)




## candidates <- stats[between(npos,100,200) & between(fsr,.25, .4) & width > 300]
## setkey(candidates,match)

## m <- candidates[2,(match)]
## reg <- get_candidate_region(m,candidates)
## reg <- GRanges(seqnames = "chr12",ranges = IRanges(start = 35968000,width = 500))

## DT <- cover_DT(reg,rep1_reads,perm = FALSE)

## M <- max(DT[,max(tags)],-DT[,min(tags)])
## M <- floor(1.2 * M)

## p <- ggplot(DT,aes(coord,tags,colour = strand))+
##   geom_step()+theme_bw()+
##   theme(legend.position = "top",plot.title = element_text(hjust  = 0))+
##     scale_color_brewer(palette = "Set1")+
##   xlab("Genomic position")+ylab("ChIP read counts")+ylim(-M,M)

## cc <- local_strand_cross_corr(rep1_reads,reg,shift = 1:150,perm = FALSE)

## s <- ggplot(cc,aes(shift,cross.corr))+geom_point(shape = 1)+
##     geom_line(size = .5,linetype = 2)+
##     geom_smooth(method = "loess",se = FALSE)+
##     geom_abline(slope = 0,intercept = 0,linetype = 2)+
##     theme_bw()+ylim(-.2,1)+
##     geom_vline(xintercept = 35,colour = "blue",linetype = 1,size = .1)+
##     geom_vline(xintercept = cc[which.max(cross.corr),(shift)],
##       colour = "red",linetype = 1,size = .1)
  
## set.seed(12345)
## seeds <- floor(runif(20) * 1e3)

## figs_dir <- "figs/local_NSC_perm"


## pdf(file = file.path(figs_dir,paste0("example_region.pdf")),height = 8, width = 6)
## grid.arrange(p + ggtitle("original"),s,nrow = 2)
## for(seed in seeds){ 
##   set.seed(seed)  
##   DU <- cover_DT(reg,rep1_reads,perm = TRUE)
##   set.seed(seed)
##   ccU <- local_strand_cross_corr(rep1_reads,reg,shift = 1:150,perm = TRUE)
##   grid.arrange(p %+% DU + ggtitle(paste("seed",seed)),s %+% ccU,nrow = 2)
## }
## dev.off()



mm <- "chr11:29512960_29513671"
reg <- get_candidate_region(mm,stats)

cc <- local_strand_cross_corr(rep1_reads,reg,shift = 1:150)

M <- 1e3
set.seed(12345)

ccs <- mclapply(1:M,function(i)
  local_strand_cross_corr(rep1_reads,reg,shift = 1:150,perm = TRUE),
  mc.cores = 16)

ccs <- mapply(function(x,y)x[,perm := paste0("M",y)],ccs,1:M,SIMPLIFY= FALSE)
ccs <- do.call(rbind,ccs)

obs <- local_strand_cross_corr(rep1_reads,reg,shift = 1:150)

dir.create(file.path(figs_dir,mm))
figs_dir <- file.path(figs_dir,mm)

## for this part we are gonna calculate summary statistics and
## calculate empirical p.values

eval_at_point <- function(s,shift,cross.corr)cross.corr[shift == s]

ss <- 34
#pdf(file = file.path(figs_dir,paste0("eval_at_",ss,".pdf")))
ggplot(ccs[,eval_at_point(ss,shift,cross.corr),by = perm],aes(V1))+
  geom_histogram(aes(y = ..density..),bins = 50,fill = NA,colour = "black")+
  stat_density(geom = "line",colour = "red")+
  geom_vline(xintercept = obs[,eval_at_point(ss,shift,cross.corr)],
    colour = "blue")
#dev.off()

## std. dev. of noise

pdf(file = file.path(figs_dir,"cross_corr_sd.pdf"))
ggplot(ccs[,sd(cross.corr),by = perm],aes(V1))+
  geom_histogram(aes(y = ..density..),bins = 50,fill = NA,colour = "black")+
  stat_density(geom = "line",colour = "red")+
  geom_vline(xintercept = obs[,sd(cross.corr)],
    colour = "blue")
dev.off()


ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

correct <- function(cover,bw)
{
  cover <- ma(as.vector(cover),n = bw)
  cover[is.na(cover)] <- 0
  return(Rle(cover))
}
 
  


grid.arrange(
s %+% local_strand_cross_corr_ma(rep1_reads,reg,shift = 1:150,bw = 0,perm = TRUE) + ggtitle("orig"),
s %+% local_strand_cross_corr_ma(rep1_reads,reg,shift = 1:150,bw = 3,perm = TRUE) + ggtitle(3),
s %+% local_strand_cross_corr_ma(rep1_reads,reg,shift = 1:150,bw = 5,perm = TRUE) + ggtitle(5),
s %+% local_strand_cross_corr_ma(rep1_reads,reg,shift = 1:150,bw = 10,perm = TRUE) + ggtitle(10),                   
s %+% local_strand_cross_corr_ma(rep1_reads,reg,shift = 1:150,bw = 15,perm = TRUE) + ggtitle(15),
nrow = 1)


local_strand_cross_corr_ma <- function (reads, region, shift = 1:300, perm = FALSE,bw = 0) 
{

  stopifnot(length(region) == 1)
  stopifnot(class(region) == "GRanges")
  stopifnot(bw >= 0)
  chr <- as.character(seqnames(region))
  stopifnot(chr %in% names(readsF(reads)) & chr %in% names(readsR(reads)))
  rF <- copy(dt2gr(readsF(reads)[[chr]]))
  rR <- copy(dt2gr(readsR(reads)[[chr]]))
  ovF <- findOverlaps(region, rF)
  ovR <- findOverlaps(region, rR)
  rF <- rF[subjectHits(ovF)]
  rR <- rR[subjectHits(ovR)]
  if (perm) {
    allr <- c(rF, rR)
    strand(allr) <- sample(as.character(strand(allr)))
    allr <- split(allr, as.character(strand(allr)))
    rF <- allr[["+"]]
    rR <- allr[["-"]]
    rm(allr)
  }
  if (length(ovF) == 0 | length(ovR) == 0) {
    cross.corr <- 0
    shift1 <- NULL
  }
  else {
    end(rF) <- start(rF)
    start(rR) <- end(rR)
    rangeF <- IRanges(min(start(rF)), max(start(rF)))
    rangeR <- IRanges(min(start(rR)), max(start(rR)))
    reg <- ranges(region)
    cF <- coverage(ranges(rF))[rangeF]
    cR <- coverage(ranges(rR))[rangeR]
    if (start(rangeF) != start(rangeR)) {
      if (start(rangeF) < start(rangeR)) {
        ext <- start(rangeR) - start(rangeF)
        cR <- c(Rle(rep(0, ext)), cR)
      }
      else {
        ext <- start(rangeF) - start(rangeR)
        cF <- c(Rle(rep(0, ext)), cF)
      }
    }
    if (end(rangeF) != end(rangeR)) {
      if (end(rangeF) < end(rangeR)) {
        ext <- end(rangeR) - end(rangeF)
        cF <- c(cF, Rle(rep(0, ext)))
      }
      else {
        ext <- end(rangeF) - end(rangeR)
        cR <- c(cR, Rle(rep(0, ext)))
      }
    }
    maxShift <- max(shift)
    if (maxShift >= length(cF)) {
      shift1 <- shift[shift < length(cF)]
    }
    else {
      shift1 <- shift
    }
    if(bw > 0){
      cF <- correct(cF,bw)
      cR <- correct(cR,bw)
    }

    if (length(shift1) > 0) {
      cc <- shiftApply(shift1, cF, cR, cor, verbose = FALSE)
    }
  }
  dt <- data.table(shift, cross.corr = 0)
  if (length(shift1) > 0) {
    dt[shift %in% shift1, `:=`(cross.corr, cc)]
  }
  dt[is.nan(cross.corr), `:=`(cross.corr, 0)]
  return(copy(dt))
}


