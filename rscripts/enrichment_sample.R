
rm(list = ls())

## the idea of this script is to compare the enriched areas across three samples

library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)

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

## get common regions to all replicates
common <- reduce(replicates[[1]]$regions)
for(k in 2:length(folders)){
  ov <- findOverlaps(common,replicates[[k]]$regions)
  common <- mapply(function(peaks,ov)peaks[queryHits(ov)],common,ov,SIMPLIFY = FALSE)
  common <- reduce(IRangesList(common))
}



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

dt2ir <- function(x)IRanges(start = x[,(start)],end = x[,(end)])


dt2gr <- function(x)
  GRanges(seqnames = x[,(seqnames)],
          ranges = dt2ir(x),
          strand = "*")

common <- as(common,"GRanges")

### subset replicates accross common peaks

replicates <- lapply(replicates,function(x,common){
  ov <- findOverlaps(dt2gr(x$stats),common)
  x[["stats"]] <- x[["stats"]][queryHits(ov)]
  setkey(x[["reads"]],match)
  common_match <- x[["stats"]][,(match)]
  x[["reads"]] <- x[["reads"]][common_match]
  return(x)},common)

library(gridExtra)
library(RColorBrewer)
library(hexbin)

enrich <- function(data,main)
{
 
  rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r <- rf(16)
  
  p1 <- ggplot(data,aes(npos,depth))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    ylim(0,3e3)+xlim(0,500)+xlab("number of unique 5' positions per region")+
    ylab("nr of reads per region")

  p2 <- ggplot(data,aes(dw_ratio,pbc))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,6)+xlab("average read coverage")+ylab("unique read coverage rate")+ylim(0,1)
  
  p3 <- ggplot(data,aes(npos,width))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,500)+ylim(0,1250)+xlab("number of unique 5' position per region")+
    ylab("width of region")

    

  p4 <- ggplot(data,aes(depth,width))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,3e3)+ylim(0,1250)+ylab("width of region")+
    xlab("nr of reads per region")

  grid.arrange(p1,p2,p3,p4,nrow = 2, ncol = 2,top = main)
  
}


gr2dt <- function(x)data.table(seqnames = as.character(seqnames(x)),start = start(x),end = end(x))
common <- gr2dt(common)
common[,match := paste(seqnames,1:length(start),sep = "_"),by = seqnames]

replicates <- lapply(replicates,function(x,common){  
  ov <- findOverlaps(dt2gr(common),dt2gr(x$stats))
  from <- x$stats[subjectHits(ov),(match)]
  to <- common[queryHits(ov),(match)]
  x$stats[,rmatch := match]
  x$stats[,match := plyr::mapvalues(match,from = from , to = to)]
  return(x)},common)
  


pdf(file = file.path(figs_dir,"enrichment_common_peaks.pdf"))
enrich(replicates[[1]]$stats,"FoxA1-rep1")
enrich(replicates[[2]]$stats,"FoxA1-rep2")
enrich(replicates[[3]]$stats,"FoxA1-rep3")
dev.off()


normalize.tagcounts <- function(counts,depth)return(counts*1e6 / depth)

Rle2dt <- function(rle_data)
{
  coord = cumsum(runLength(rle_data))
  counts = runValue(rle_data)
  return(data.table(coord=coord,counts = counts))
}



setkey(common,match)
setkey(replicates[[1]]$stats,match)
setkey(replicates[[2]]$stats,match)
setkey(replicates[[3]]$stats,match)
setkey(replicates[[1]]$reads,match)
setkey(replicates[[2]]$reads,match)
setkey(replicates[[3]]$reads,match)



cover_dt <- function(reads,key,depth,st,fl)
{
  cover <- coverage(resize(dt2ir(reads[key][strand == st]),fl))
  out <- Rle2dt(cover)
  out[,strand := st]
  if(st == "+"){
    sgn <- 1
  }else{
    sgn <- -1
  }
  out[,counts := sgn * normalize.tagcounts(counts,depth)]
  return(out)
}

fill_cover_dt <- function(repl,covers,base)
{
  cover <- covers[[repl]]
  base2 <- copy(base)
  if(nrow(cover) > 0){
    base2[ between(coord,cover[,min(coord)],cover[,max(coord)]),
      counts := approxfun(x = cover[,(coord)],y = cover[,(counts)],
        method = "constant")(coord)]
  }
  base2[,rep := repl]
  return(base2)
}

get_profile <- function(key,common,replicates,fl = 1,depth = TRUE)
{
    
  chr <- common[key,(seqnames)]
  start <- common[key,(start)]
  end <- common[key,(end)]

  rkeys <- lapply(replicates,function(x,key)x$stats[key,(rmatch)],key)
  reads <- lapply(replicates,function(x)x$reads)

  if(!depth){
    depths <- rep(1,length(replicates))
  }else{
    depths <- sapply(replicates,function(x)x$depth)

  }

 
  fwds <- mapply(cover_dt,reads,rkeys,depths,
    MoreArgs = list(st = "+",fl = fl),SIMPLIFY = FALSE)
  bwds <- mapply(cover_dt,reads,rkeys,depths,
   MoreArgs = list(st = "-",fl = fl),SIMPLIFY = FALSE)

  fwds2 <- data.table(coord = start:end,counts = 0 ,strand = "+")
  bwds2 <- data.table(coord = start:end,counts = 0 ,strand = "-")

  fwds2 <- lapply(1:length(replicates),fill_cover_dt,fwds,fwds2)
  fwds2 <- do.call(rbind,fwds2)

  bwds2 <- lapply(1:length(replicates),fill_cover_dt,bwds,bwds2)
  bwds2 <- do.call(rbind,bwds2)

  M <- max(fwds2[,(counts)],-bwds2[,(counts)])

  p = ggplot(data.frame(x=0,y=0))+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = I("black"))+
    geom_step(data = fwds2,aes(x=coord,y=counts),colour = "red",direction = "vh")+
    geom_step(data = bwds2,aes(x=coord,y=counts),colour = "blue",direction = "vh")+
    facet_grid( rep ~ . ,scales  = "free_y")+
    theme(legend.position = "none")+
    ylab("Normalized counts")+
    scale_x_continuous(limits = c(start , end))+
    ggtitle(paste(chr,":",prettyNum(start,big.mark = ","),"-",prettyNum(end,big.mark = ",")))

## scale_y_continuous( limits = 1.2*M * c(-1,1))+
  
  return(p)
}



## we think that the best quality set comes from rep1
## therefore we are going to look by their characteristics

## dw_ratio > 2

plot_keys <- function(keys,common,replicates)
{
  for(k in 1:length(keys)){
    message(k)
    p <- get_profile(keys[k],common,replicates)
    print(p)
  }
}

stats <- lapply(replicates,function(x)x$stats)

reads <- lapply(replicates,function(x)x$reads)
setkey(reads[[1]],match)
setkey(reads[[2]],match)
setkey(reads[[3]],match)




keys <- replicates[[1]]$stats[dw_ratio > 2,(match)]
pdf(file = file.path(figs_dir,"dw_ratio_gr_2.pdf"),width = 6,height = 9)
plot_keys(keys,common,replicates)
dev.off()




## all samples dw_ratio
keys <- lapply(stats,function(x)x[dw_ratio > 2,(match)])
keys <- intersect(keys[[1]],intersect(keys[[2]],keys[[3]]))
pdf(file = file.path(figs_dir,"dw_ratio_gr_2_all.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


set.seed(123321)

## rep1, npos > 200
keys <- stats[[1]][ npos > 200,sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_gr_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()



## rep2, npos > 200
keys <- stats[[2]][ npos > 200,(match)]
pdf(file = file.path(figs_dir,"rep2_npos_gr_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


## rep3, npos > 200
keys <- stats[[3]][ npos > 200,(match)]
pdf(file = file.path(figs_dir,"rep3_npos_gr_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


set.seed(123321)

## rep1, npos between 100 and 200
keys <- stats[[1]][ between(npos,100,200),sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_bet_100_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

## rep1, npos between 100 and 200
keys <- stats[[2]][ between(npos,100,200),sample(match,500)]
pdf(file = file.path(figs_dir,"rep2_npos_bet_100_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

## rep3 npos between 100 and 200
keys <- stats[[3]][ between(npos,100,200),(match)]
pdf(file = file.path(figs_dir,"rep3_npos_bet_100_200.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

set.seed(123321)

## rep1, npos between 50 and 100
keys <- stats[[1]][ between(npos,50,100),sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_bet_50_100.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


keys <- stats[[2]][ between(npos,50,100),sample(match,500)]
pdf(file = file.path(figs_dir,"rep2_npos_bet_50_100.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

keys <- stats[[3]][ between(npos,50,100),(match)]
pdf(file = file.path(figs_dir,"rep3_npos_bet_50_100.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

set.seed(123321)

## rep1, npos between 25 and 50
keys <- stats[[1]][ between(npos,25,50),sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_bet_25_50.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


## rep2, npos between 25 and 50
keys <- stats[[2]][ between(npos,25,50),sample(match,500)]
pdf(file = file.path(figs_dir,"rep2_npos_bet_25_50.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

## rep3, npos between 25 and 50
keys <- stats[[3]][ between(npos,25,50),sample(match,500)]
pdf(file = file.path(figs_dir,"rep3_npos_bet_25_50.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


set.seed(123321)

## rep1, npos between 10 and 25
keys <- stats[[1]][ between(npos,10,25),sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_bet_10_25.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


## rep2, npos between 10 and 25
keys <- stats[[2]][ between(npos,10,25),sample(match,500)]
pdf(file = file.path(figs_dir,"rep2_npos_bet_10_25.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)

dev.off()

## rep3, npos between 10 and 25
keys <- stats[[3]][ between(npos,10,25),sample(match,500)]
pdf(file = file.path(figs_dir,"rep3_npos_bet_10_25.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

set.seed(123321)

## rep1, npos between 5 and 10
keys <- stats[[1]][ between(npos,5,10),sample(match,500)]
pdf(file = file.path(figs_dir,"rep1_npos_bet_5_10.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

## rep1, npos between 5 and 10
keys <- stats[[2]][ between(npos,5,10),sample(match,500)]
pdf(file = file.path(figs_dir,"rep2_npos_bet_5_10.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()

## rep1, npos between 5 and 10
keys <- stats[[3]][ between(npos,5,10),sample(match,500)]
pdf(file = file.path(figs_dir,"rep3_npos_bet_5_10.pdf"),width = 6 ,height = 9)
plot_keys(keys,common,replicates)
dev.off()


### try to build a signal - to - noise measure
stats <- lapply(replicates,function(x)x$stats)

reads <- lapply(replicates,function(x)x$reads)
setkey(reads[[1]],match)
setkey(reads[[2]],match)
setkey(reads[[3]],match)

keys <- stats[[1]][dw_ratio > 2,(match)]
kk <- sample(keys,1)

### get reads


snr <- function(key,stats,reads,repl = 1)
{
  browser()
  repl_key <- stats[[repl]][key,(rmatch)]
  peak_reads <- copy(reads[[repl]][repl_key])
  setkey(peak_reads,strand)
  fwd <- peak_reads["+"]
  bwd <- peak_reads["-"]
  

}


par(mfrow = c(2,1))
plot(table(fwd[,sort(start)]))
##plot(acf(table(fwd[,sort(start)])))
plot(table(bwd[,sort(end)]))
##plot(acf(table(bwd[,sort(end)])))
dev.off()

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}



snr(kk,stats,reads,1)

stats <- lapply(replicates,function(x)x$stats)

reads <- lapply(replicates,function(x)x$reads)
setkey(reads[[1]],match)
setkey(reads[[2]],match)
setkey(reads[[3]],match)




keys <- stats[[1]][seqnames != "chrM" & npos > 100 ,(match)]
kk <- sample(keys,1)

pdf(file = "example2.pdf")
reads1 <- copy(reads[[1]][ stats[[1]][kk,(rmatch)]])
reads1[,end := end - min(start)]
reads1[,start := start - min(start)]
reads2 <- copy(reads[[2]][ stats[[2]][kk,(rmatch)]])
reads2[,end := end - min(start)]
reads2[,start := start - min(start)]
reads3 <- copy(reads[[3]][ stats[[3]][kk,(rmatch)]])
reads3[,end := end - min(start)]
reads3[,start := start - min(start)]
get_profile(kk,common,replicates,depth = TRUE)
cc1 <- as.numeric(coverage(dt2ir(reads1)))
cc2 <- as.numeric(coverage(dt2ir(reads2)))
cc3 <- as.numeric(coverage(dt2ir(reads3)))
plot(cc1,type = "l")
hist(cc1,breaks = 100)
plot(cc2,type = "l")
hist(cc2,breaks = 100)
plot(cc3,type = "l")
hist(cc3,breaks = 100)
dev.off()






