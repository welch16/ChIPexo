
rm(list = ls())

## the idea of this script is to compare the enriched areas across three samples

library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)

folders <- paste0("FoxA1-rep",1:3)
dir <- "/p/keles/ChIPexo/volume3/Analysis/Carroll/mouse"
figs_dir <- "figs/local_cross_corr/all_reps"


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


## we think that the best quality set comes from rep1
## therefore we are going to look by their characteristics

## dw_ratio > 2

stats <- lapply(replicates,function(x)x$stats)

reads <- lapply(replicates,function(x)x$reads)
setkey(reads[[1]],match)
setkey(reads[[2]],match)
setkey(reads[[3]],match)


library(devtools)
load_all("~/Desktop/Docs/Code/ChIPUtils")

N <- 300
set.seed(123321)
keys1 <- stats[[1]][ npos > 200,sample(match,N)]
set.seed(123321)
keys2 <- stats[[1]][ between(npos,100,200),sample(match,N)]
set.seed(123321)
keys3 <- stats[[1]][ between(npos,50,100),sample(match,N)]
set.seed(123321)
keys4 <- stats[[1]][ between(npos,25,50),sample(match,N)]

create_reads <- function(rr)
{
  if(nrow(rr) > 0){    
    setkey(rr,strand)
    if(nrow(rr["+",nomatch = 0]) == 0 | nrow(rr["-",nomatch = 0]) == 0){
      out <- NA
    }else{
      out <- new("reads",readsFile = "",readsF = split(rr["+"],rr["+",(seqnames)]),
      readsR = split(rr["-"],rr["-",(seqnames)]),nReads = nrow(rr))
    }
  }else{
    out <- NA
  }
  return(out)
}


build_reads <- function(key,common,reads)
{

  rkeys <- lapply(replicates,function(x,key)x$stats[key,(rmatch)],key)
  reads <- lapply(replicates,function(x)x$reads)

  for(k in 1:length(reads)){
    setkey(reads[[k]],match)
  }
 
  subreads <- mapply(function(x,r)x[r,nomatch = 0],reads,rkeys,SIMPLIFY = FALSE)
  subreads <- lapply(subreads,create_reads)

  return(subreads)

}

reads1 <- mclapply(keys1,build_reads,common,reads,mc.cores = mc)
reads2 <- mclapply(keys2,build_reads,common,reads,mc.cores = mc)
reads3 <- mclapply(keys3,build_reads,common,reads,mc.cores = mc)

get_local_scc <- function(key,rr,common,shift)
{

  reg <- dt2gr(common[key])

  idx <- !is.na(rr)
  cross_corr <- lapply(rr[idx],local_strand_cross_corr,reg,shift)

  if( any(!idx)){
    cross_corr[[names(which(!idx))]] <- data.table(shift,cross.corr = NA)
  }
  
  cross_corr <- mapply(function(x,y)x[,sample:= y],cross_corr,names(cross_corr),
    SIMPLIFY = FALSE)
  cross_corr <- do.call(rbind,cross_corr)
  
  return(cross_corr)

}

shift <- 1:200
cc1 <- mcmapply(get_local_scc,keys1,reads1,MoreArgs = list(common,shift),SIMPLIFY = FALSE,
         mc.cores = mc,mc.preschedule = TRUE)

cc2 <- mcmapply(get_local_scc,keys2,reads2,MoreArgs = list(common,shift),SIMPLIFY = FALSE,
         mc.cores = mc,mc.preschedule = TRUE)
cc3 <- mcmapply(get_local_scc,keys3,reads3,MoreArgs = list(common,shift),SIMPLIFY = FALSE,
         mc.cores = mc,mc.preschedule = TRUE)

get_plot <- function(cc)
{
  p <- ggplot(cc,aes(shift,cross.corr))+geom_point(size = 1.3)+
    geom_line(linetype = 3,size = .25)+
    geom_smooth(method = "loess",se = FALSE)+
    facet_grid(sample ~ . )
  return(p)
}

plots1 <- mclapply(cc1,get_plot,mc.cores = mc)
plots2 <- mclapply(cc2,get_plot,mc.cores = mc)
plots3 <- mclapply(cc3,get_plot,mc.cores = mc)

## rep1, npos > 200
pdf(file = file.path(figs_dir,"rep1_npos_gr_200_cc.pdf"),width = 6 ,height = 9)
u <- lapply(plots1,print)
dev.off()


## rep1, npos between 100 and 200
pdf(file = file.path(figs_dir,"rep1_npos_bet_100_200_cc.pdf"),width = 6 ,height = 9)
u <- lapply(plots2,print)
dev.off()



## rep1, npos between 50 and 100
pdf(file = file.path(figs_dir,"rep1_npos_bet_50_100_cc.pdf"),width = 6 ,height = 9)
u <- lapply(plots3,print)
dev.off()

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



s1 <- lapply(cc1,function(x)x[,noise(shift,cross.corr),by = sample])
s2 <- lapply(cc2,function(x)x[,noise(shift,cross.corr),by = sample])
s3 <- lapply(cc3,function(x)x[,noise(shift,cross.corr),by = sample])
  

s1 <- mapply(function(x,y){
  x[,name := y]
  return(x)},s1,names(s1),SIMPLIFY = FALSE)
s1 <- do.call(rbind,s1)

s2 <- mapply(function(x,y){
  x[,name := y]
  return(x)},s2,names(s2),SIMPLIFY = FALSE)
s2 <- do.call(rbind,s2)

s3 <- mapply(function(x,y){
  x[,name := y]
  return(x)},s3,names(s3),SIMPLIFY = FALSE)
s3 <- do.call(rbind,s3)

s1[,block := "npos > 200"]
s2[,block := "100 < npos < 200"] 
s3[,block := "50 < npos < 100"] 

ss <- rbind(s1,s2,s3)
ss[ , block := factor(block ,levels = rev(c("50 < npos < 100","100 < npos < 200","npos > 200")))]

pdf(file = file.path(figs_dir,"SCC_loc_poly_reg_noise.pdf"))
ggplot(ss , aes(block, V1,colour = sample))+geom_boxplot()+
  facet_grid(. ~ sample  )+ylim(0.02,.1)+theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))+
  scale_color_brewer(palette = "Set1")+ylab("noise")+xlab("")
dev.off()

## can we approximate nsc as cc[fl] / noise ?

nsc1 <- lapply(cc1,function(x)x[,max(cross.corr),by = sample])
nsc2 <- lapply(cc2,function(x)x[,max(cross.corr),by = sample])
nsc3 <- lapply(cc3,function(x)x[,max(cross.corr),by = sample])


nsc1 <- mapply(function(x,y){
  x[,name := y]
  return(x)},nsc1,names(nsc1),SIMPLIFY = FALSE)
nsc1 <- do.call(rbind,nsc1)

nsc2 <- mapply(function(x,y){
  x[,name := y]
  return(x)},nsc2,names(nsc2),SIMPLIFY = FALSE)
nsc2 <- do.call(rbind,nsc2)

nsc3 <- mapply(function(x,y){
  x[,name := y]
  return(x)},nsc3,names(nsc3),SIMPLIFY = FALSE)
nsc3 <- do.call(rbind,nsc3)

nsc1[,block := "npos > 200"]
nsc2[,block := "100 < npos < 200"] 
nsc3[,block := "50 < npos < 100"] 

nsc <- rbind(nsc1,nsc2,nsc3)

pdf(file = file.path(figs_dir,"argmax_scc.pdf"))
ggplot(nsc , aes(block, V1,colour = sample))+geom_boxplot()+
  facet_grid(. ~ sample  )+theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))+
  scale_color_brewer(palette = "Set1")+ylab("max scc")+xlab("")
dev.off()

nsc <- nsc[, s := ss[,(V1)]]
nsc <- nsc[, NSC := V1 / sqrt(s)]

pdf(file = file.path(figs_dir,"argmax_scc.pdf"))

ggplot(nsc , aes(block, NSC,colour = sample))+geom_boxplot()+
  facet_grid(. ~ sample  )+theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))+ylim(0,5)+
  scale_color_brewer(palette = "Set1")+ylab("NSC")+xlab("")
dev.off()
