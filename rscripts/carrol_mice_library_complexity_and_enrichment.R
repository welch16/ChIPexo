
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

pbc <- sapply(reads,PBC)
depth <- sapply(reads,nreads)
ssd <- sapply(reads,SSD) ## there is a small issue with this function. sometimes can get numeric overflow
## Warning messages:
## 1: In sum(w) : integer overflow - use sum(as.numeric(.))
## 2: In sum(w) : integer overflow - use sum(as.numeric(.))
## 3: In sum(w) : integer overflow - use sum(as.numeric(.))



FSR <- function(reads)
{
  summ <- summary(reads)
  rF <- summ[,sum(readsF)]
  rR <- summ[,sum(readsR)]

  return(rF / (rF + rR))
}

fsr <- sapply(reads,FSR)

exo <- data.table(sample  = names(reads))

exo[,nreads := depth]
exo[,pbc := pbc]
exo[,ssd := ssd] 
exo[,fsr := fsr]

mouse.size <- data.table(read.table(
 system.file( "extdata","chrom.sizes","mm9.chrom.sizes",package = "ChIPUtils")))

scc <- lapply(reads,strand_cross_corr,shift = 1:300,chrom.sizes = mouse.size,parallel = TRUE)

nsc <- sapply(scc,function(x)x[ ,max(cross.corr) / min(cross.corr)])
exo[, nsc := nsc]

## Carroll mouse
## ERR336942.bam FoxA1-rep1 
## ERR336956.bam FoxA1-rep2
## ERR336935.bam FoxA1-rep3

bamfiles <-c("ERR336942.bam","ERR336956.bam","ERR336935.bam")
repl <- paste0("rep-",1:3)

scc <- mapply(function(x,y)x[,sample := y],scc,exo[,(sample)],SIMPLIFY = FALSE)
SCC <- do.call(rbind,scc)
SCC[ , sample := plyr::mapvalues(sample ,from = bamfiles ,to =repl)]

exo[,which.max := SCC[,which.max(cross.corr),by = sample][,(V1)]]

exo[, sample := plyr::mapvalues(sample ,from = bamfiles ,to =repl)]

save(exo,file = "data/for_paper/Carroll_mouse_FoxA1_summary.RData")

figs_dir <- "figs/Carroll_mice_for_paper/"

pdf(file = file.path(figs_dir,"Strand_cross_corr.pdf"),width = 9 , height = 5)
ggplot(SCC,aes(shift,cross.corr,colour = sample))+geom_line(size = .8)+
  scale_color_brewer(name = "",palette = "Dark2")+theme_bw()+theme(legend.position = "top")+
  ylab("Strand cross-correlation")
dev.off()


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

common <- subsetByOverlaps(dt2gr(regs[[1]]),subsetByOverlaps(dt2gr(regs[[2]]),dt2gr(regs[[3]])))

build_stats <- function(region,reads)
{
  ## fix formats and stuff
  
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

  stats[ , M := as.numeric(NA)]
  stats[ , A := as.numeric(NA)]

  stats[f > 0 & r > 0, M := log2(f * r) - 2 * log2(width)]
  stats[f > 0 & r > 0, A := log2( f/ r)]

  stats[ , strand := NULL]

  return(stats)
}

regs <- lapply(regs,dt2gr)

stats <- mcmapply(build_stats,regs,gr,mc.cores = 6 ,SIMPLIFY = FALSE)

stats_common <- mclapply(gr,function(x)build_stats(common,x),mc.cores = 6)

aux <- mapply(function(x,y)x[,sample := y],stats,names(stats),SIMPLIFY =FALSE)
aux <- do.call(rbind,aux)

pdf(file = file.path(figs_dir,"FoxA1_number_unique_positions_perSample.pdf"))
ggplot(aux[npos > 10] , aes(npos))+  geom_histogram()+scale_x_log10()+xlim(0,100)+facet_wrap( ~  sample)
dev.off()

library(hexbin)
library(scales)
r <- viridis::viridis(100,option = "D")

pdf(file = file.path(figs_dir,"FoxA1_enrichment.pdf"),width = 9 , height = 5 )
p <- ggplot(aux , aes( ave_reads,cover_rate))+stat_binhex(bins = 50)+
  facet_wrap( ~ sample)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,4)+ylim(0,1)+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Average read coverage")+
  ylab("Unique read coverage rate")
print(p + ggtitle("A"))
print(p + ggtitle("All regions"))
print( p  %+% aux[ npos > 10] + ggtitle("A"))  ## npos > 10 
print( p  %+% aux[ npos > 30] + ggtitle("B")) ## npos > 30
dev.off()

enrichment <- list( data = aux , plot = p)
save(enrichment , file = "data/for_paper/FoxA1_enrich_plot.RData")


## signal-to-noise

### filter to regions where we can calcualte the local scc

pdf(file = file.path(figs_dir,"FoxA1_MA_plots.pdf"),width = 9 , height = 5 )
p <-ggplot(aux[f > 0 & r > 0 ] , aes( M , A))+stat_binhex(bins = 50)+
  facet_wrap( ~ sample)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))
print(p + ggtitle("All regions"))
print( p  %+% aux[ npos > 10] + ggtitle("Npos > 10")) 
print( p  %+% aux[ npos > 30] + ggtitle("Npos > 30")) 
dev.off()


filter_stats <- lapply(stats,function(x)x[ f > 0 & r > 0])

local_scc <- function(stat , reads , min_npos,max_npos, nsample, shift = 1:300)
{

  stat <- copy(stat[ between(npos, min_npos,max_npos)])
  message(nrow(stat))
  if(nrow(stat) > nsample){
    stat <- stat[sample(match , nsample)]
  }

  regions <- dt2gr(stat[,2:4,with = FALSE])
  regions <- split(regions,start(regions))
  names(regions) <- stat[,(match)]

  regions <- as.list(regions)
  
  scc <- mclapply(regions,function(x){
    out <- local_strand_cross_corr(reads,x,shift)
    return(out)},mc.cores = 20)

  return(scc)
}

nsamp <- 400
shift <- 1:200

strata_high <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 100, max_npos = Inf,nsample = nsamp ,shift = shift),SIMPLIFY = FALSE)

strata_med <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 50, max_npos = 100,nsample = nsamp ,shift = shift),SIMPLIFY = FALSE)

strata_low <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 20, max_npos = 50,nsample = nsamp ,shift = shift),SIMPLIFY = FALSE)


local_scc_data <- list(strata_high, strata_med,strata_low)

save(local_scc_data , file = "data/for_paper/Carroll_FoxA1_local_scc_by_strata.RData")

## get it as a big table
strati <- c("high","med","low")
load(file = "data/for_paper/Carroll_FoxA1_local_scc_by_strata.RData")


local_scc_data <- lapply(local_scc_data,function(x,nm){
  names(x) <- nm
  return(x)},bamfiles)

names(local_scc_data) <- strati

local_scc_data <- lapply(local_scc_data,
  function(x){
    out <- lapply(x,function(y){
    scc <- mapply(function(curve,name)curve[,match := name],y,names(y),SIMPLIFY = FALSE)
    scc <- do.call(rbind,scc)
    return(scc)
    })
    out <- mapply(function(scc,name)scc[,file := name],out,names(out),SIMPLIFY = FALSE)
    out <- do.call(rbind,out)
    return(out)
  })

local_scc_data <- mapply(function(x,y)x[,strata := y],local_scc_data,strati,SIMPLIFY = FALSE)
local_scc_data <- do.call(rbind,local_scc_data)

noise <- function(shift,cross.corr)
{
  if(all( is.na(cross.corr))){
    out <- Inf
  }else{
    mod <- loess(cross.corr ~ shift)
    out <- mod$s
  }
  return(out)
}

cc_max <- function(shift,cross.corr)
{
  if(all(is.na(cross.corr))){
    out <- Inf
  }else{
    
    mod <- loess(cross.corr ~ shift)
    out <- max(predict(mod))
  }
  return(out)
}
 

ss <- local_scc_data[,noise(shift,cross.corr), by = .(match, file,strata)]
nsc <- local_scc_data[,cc_max(shift,cross.corr), by = .(match, file, strata)]

setnames(ss , names(ss), c("match","file","strata","noise"))
setnames(nsc , names(nsc), c("match","file","strata","max"))

nsc[ , file := NULL]
nsc[ , strata := NULL]

dat <- merge(ss,nsc,by = "match",allow.cartesian = TRUE)

dat[ , nsc := max / noise ]

dat[ , strata := factor(strata, levels = rev(strati) ) ]

## dat[ , strata := plyr::mapvalues(strata ,
##          from = c(
##            "high",
##            "med1",
##            "med2",
##            "low1",
##            "low2",
##            "low3"),
##          to = c(
##            "(100,Inf)",
##            "(75,100)",
##            "(50,75)",
##            "(30,50)",
##            "(15,30)",
##            "(10,15)"             
##            ))]

strata_plots <- list()
strata_plots[[1]] <- ggplot(dat , aes( strata , noise,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.25)
strata_plots[[2]] <- ggplot(dat , aes( strata , max,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Loess max")+geom_abline( slope = 0 , intercept = 0, linetype = 2)+
  ylim(0,.3)
strata_plots[[3]] <- ggplot(dat , aes( strata , nsc,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(0,3)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("B")
strata_plots[[4]] <- ggplot(dat , aes( file , nsc,colour = strata))+geom_boxplot()+facet_grid( . ~ strata) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(0,3)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("B")

pdf(file = file.path(figs_dir,"Local_SCC_indicator_by_strata.pdf"),width = 9,height = 4)
u <- lapply(strata_plots,print)
dev.off()



anchor_stats <- filter_stats[[2]]

nsamp <- 400
shift <- 1:200

high <- anchor_stats[ npos > 100][ sample(match,nsamp)]
med <- anchor_stats[between(npos,50,100)][sample(match,nsamp)]
low <- anchor_stats[between(npos,20,50)][sample(match,nsamp)]

high_stats <- mclapply(gr,function(x)build_stats(dt2gr(high[,2:4,with = FALSE]),x),mc.cores =3)
med_stats <- mclapply(gr,function(x)build_stats(dt2gr(med[,2:4,with = FALSE]),x),mc.cores =3)
low_stats <- mclapply(gr,function(x)build_stats(dt2gr(low[,2:4,with = FALSE]),x),mc.cores =3)

high_stats <- lapply(high_stats,function(x)x[f > 0 & r > 0])
med_stats <- lapply(med_stats,function(x)x[f > 0 & r > 0])
low_stats <- lapply(low_stats,function(x)x[f > 0 & r > 0])

local_scc2 <- function(stat , reads , shift = 1:300)
{

  regions <- dt2gr(stat[,2:4,with = FALSE])
  regions <- split(regions,start(regions))
  names(regions) <- stat[,(match)]

  regions <- as.list(regions)
  
  scc <- mclapply(regions,function(x){
    out <- local_strand_cross_corr(reads,x,shift)
    return(out)},mc.cores = 20)

  return(scc)
}

strata_high2 <- mapply(local_scc2,high_stats,reads,
      MoreArgs = list(shift = shift),SIMPLIFY = FALSE)

strata_med2 <- mapply(local_scc2,med_stats,reads,
      MoreArgs = list(shift = shift),SIMPLIFY = FALSE)

strata_low2 <- mapply(local_scc2,low_stats,reads,
      MoreArgs = list(shift = shift),SIMPLIFY = FALSE)

local_scc_data <- list(strata_high2,strata_med2,strata_low2)
names(local_scc_data) <- c("high","med","low")

local_scc_data <- lapply(local_scc_data,
  function(x){
    out <- lapply(x,function(y){
    scc <- mapply(function(curve,name)curve[,match := name],y,names(y),SIMPLIFY = FALSE)
    scc <- do.call(rbind,scc)
    return(scc)
    })
    out <- mapply(function(scc,name)scc[,file := name],out,names(out),SIMPLIFY = FALSE)
    out <- do.call(rbind,out)
    return(out)
  })

local_scc_data <- mapply(function(x,y)x[,strata := y],local_scc_data,strati,SIMPLIFY = FALSE)
local_scc_data <- do.call(rbind,local_scc_data)



ss <- local_scc_data[,noise(shift,cross.corr), by = .(match, file,strata)]
nsc <- local_scc_data[,cc_max(shift,cross.corr), by = .(match, file, strata)]

setnames(ss , names(ss), c("match","file","strata","noise"))
setnames(nsc , names(nsc), c("match","file","strata","max"))

nsc[ , file := NULL]
nsc[ , strata := NULL]
dat2 <- merge(ss,nsc,by = "match",allow.cartesian = TRUE)
dat2[ , nsc := max / noise ]
dat2[ , strata := factor(strata, levels = rev(strati) ) ]

dat2[  , file := gsub(".bam","",file)]

strata_plots <- list()
strata_plots[[1]] <- ggplot(dat2 , aes( strata , noise,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.25)
strata_plots[[2]] <- ggplot(dat2 , aes( strata , max,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Loess max")+geom_abline( slope = 0 , intercept = 0, linetype = 2)+
  ylim(0,.3)
strata_plots[[3]] <- ggplot(dat2 , aes( strata , nsc,colour = file))+geom_boxplot()+facet_grid( . ~ file) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(0,3)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("B")
strata_plots[[4]] <- ggplot(dat2 , aes( file , nsc,colour = file))+geom_boxplot()+facet_grid( . ~ strata) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(0,2)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("C")

pdf(file = file.path(figs_dir,"Local_SCC_indicator_by_strata_sameRegion.pdf"),width = 9,height = 4)
u <- lapply(strata_plots,print)
dev.off()

library(grid)
library(gridExtra)

plots <- list()
plots[[1]] <- ggplot(dat , aes( file , noise , colour = file))+geom_boxplot()+
  facet_grid( . ~ file ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.18)
plots[[2]] <- ggplot(dat , aes( file , max , colour = file))+geom_boxplot()+
  facet_grid( . ~ file ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Loess max")+ylim(-.05, .4)+geom_abline(slope = 0,intercept = 0,linetype = 2)
plots[[3]] <- ggplot(dat , aes( file, nsc , colour = file))+geom_boxplot()+
  facet_grid( . ~ file ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Local NSC")+ylim(-1, 3)+geom_abline(slope = 0,intercept = 0,linetype = 2)


pdf(file = file.path(figs_dir,"Local_SCC_indicators.pdf"),width = 5 , height = 6)
u <- lapply(plots,print)
dev.off()




### example_plot

segvis_to_go <- function(reads,region)
{
  fwd <- dt2gr(readsF(reads)[[1]])
  end(fwd) <- start(fwd)
  bwd <- dt2gr(readsR(reads)[[1]])
  start(bwd) <- end(bwd)
  
  fwd <- subsetByOverlaps(fwd ,region)
  bwd <- subsetByOverlaps(bwd ,region)

  fwd <- coverage(ranges(fwd))
  bwd <- coverage(ranges(bwd))

  dt <- data.table( coord = start(region):end(region))
  dt[ , fwd := stepfun(cumsum(runLength(fwd)),c(0,runValue(fwd)))(coord)]
  dt[ , bwd := -stepfun(cumsum(runLength(bwd)),c(0,runValue(bwd)))(coord)]

  dt <- melt(dt, id.vars = "coord")
    
  return(dt)
  
}



  out <- ggplot(dt , aes(coord ,value , colour = variable ))+geom_step()+
    scale_color_brewer(palette = "Set1",name = "")+theme_bw()+
      theme(legend.position = "none")+geom_abline(slope = 0 ,intercept = 0 , linetype =2)

