
rm(list = ls())

## the idea of this script is to compare the enriched areas across three samples

library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)

folders <- paste0("FoxA1-rep",1:3)
dir <- "/p/keles/ChIPexo/volume3/Analysis/Carroll/mouse"


build_files <- function(folder,dir)
{
  files <- list.files(file.path(dir,folder,"data"))
  patterns <- c("depth","reads_by_region","regions","summary_stat")
  files <- files[sapply(patterns,function(x,files)
                        grep(x,files),files)]
  return(files)
}

files <- lapply(folders,build_files,dir)
mc <- 8


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


## tidy
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

replicates <- lapply(replicates,function(X){
  summ <- mclapply(X[["summary"]],
  function(x){
    x[,region := paste(chrID,region,sep = "_")]
    x[,chrID := NULL]
    nms <- names(x)
    setnames(x,nms,c("match",nms[-1]))
    return(x)},mc.cores = mc)
  summ <- do.call(rbind,summ)
  X[["summary"]] <- summ
  return(X)})

replicates <- lapply(replicates,function(x){
  reg <- x[["regions"]]
  summ <- x[["summary"]]
  xx <- merge(reg,summ,by = "match",all = TRUE)
  x[["regions"]] <- NULL
  x[["summary"]] <- xx
  return(x)})
                     

common <- mapply(function(x,y){
  GRanges(seqnames = y, ranges = x ,strand = "*")},
  common,names(common),SIMPLIFY = FALSE)                 
common <- GRangesList(common)

common <- unlist(common)
gr2dt <- function(x)data.table(seqnames = as.character(seqnames(x)),start = start(x),end = end(x))
common <- gr2dt(common)
common[,match := paste(seqnames,1:length(start),sep = "_"),by = seqnames]

setkey(replicates[[1]]$reads,match)
setkey(replicates[[2]]$reads,match)
setkey(replicates[[3]]$reads,match)



create_reads_from_key <- function(key,rr)
{ 
  qq <- copy(rr[key])
  qq[,match := NULL]
  setkey(qq,strand)
  return(new("reads",readsFile = "",readsF = split(qq["+"],qq["+",(seqnames)]),
    readsR = split(qq["-"],qq["-",(seqnames)]),nReads = nrow(qq)))             
}




## repl <- replicates[[1]]
## matchs <- repl$summary[label == "both" & npos > 50,(match)]

## cc = mclapply(matchs[1:2000],loc_scc_wrap,repl,shift = seq(0,100,by = 1),mc.cores = mc)

## lapply(cc[which(sapply(cc,class) != "try-error")] ,function(x){
##   plot(x)
##   lines(x)
##       })
## dev.off()


library(devtools)
load_all("~/Desktop/Docs/Code/ChIPUtils")


stats <- lapply(replicates,function(x)x[["summary"]])

loc_scc_wrap <- function(key,repl,shift = 1:100)
{
  message(key)
  reg <- repl[["summary"]][key]
  reads <- create_reads_from_key(key,repl[["reads"]])
  reg <- dt2gr(reg[,2:4,with = FALSE])
  tryCatch(local_strand_cross_corr(reads,reg,shift),
           warning = function(w) NA,
           error = function(e)NA,
           finally = NA)
}

plot_cc <- function(cc1,key,stat)
{
  if(!is.na(cc1)){
    st <- stat[key]
    tt <- paste(st[,(seqnames)],":",prettyNum(st[,(start)],big.mark = ","),"-",
                prettyNum(st[,(end)],big.mark = ","),sep = " ") 
    p <- ggplot(cc1,aes(shift, cross.corr))+geom_point(size = 1.2)+geom_line(size = .3,linetype = 2)+
      geom_smooth(method = "loess",se = FALSE)+ggtitle(tt)
    return(p)
  }
}


noise_est <- function(cc1)
{
  if(is.na(cc1)){
    out <- NA
  }else{
    summ <- cc1[,summary(loess(cross.corr ~ shift))]
    out <- summ$s
  }
  return(out)
}


figs_dir <- "figs/local_cross_corr"

## rep1, npos > 200
set.seed(123321)
keys <- stats[[1]][ npos > 200,sample(match,500)]
cc1 <- mclapply(keys,loc_scc_wrap,replicates[[1]],shift = 1:200,mc.cores = mc)
s1 <- do.call(c,mclapply(cc1,noise_est,mc.cores = mc))





pdf(file  = file.path(figs_dir,"rep1_npos_gr_200.pdf"))
u <- mcmapply(plot_cc,cc1,keys,MoreArgs = list(stats[[1]]),mc.cores = mc,SIMPLIFY = FALSE)
Z <- lapply(u,print)
dev.off()

 


## rep1, npos between 100 and 200
set.seed(123321)
keys <- stats[[1]][ between(npos,100,200),sample(match,500)]
cc2 <- mclapply(keys,loc_scc_wrap,replicates[[1]],shift = 1:200,mc.cores = mc)
s2 <- do.call(c,mclapply(cc2,noise_est,mc.cores = mc))


pdf(file = file.path(figs_dir,"rep1_npos_bet_100_200.pdf"),width = 6 ,height = 9)
u <- mcmapply(plot_cc,cc,keys,MoreArgs = list(stats[[1]]),mc.cores = mc,SIMPLIFY = FALSE)
Z <- lapply(u,print)
dev.off()



## rep1, npos between 50 and 100
set.seed(123321)
keys <- stats[[1]][ between(npos,50,100),sample(match,500)]
cc3 <- list()
for(k in keys){
  if( stats[[1]][k , (label) ]!= "both"){
    cc3[[k]] <- NA
  }else{
    cc3[[k]] <- loc_scc_wrap(k,replicates[[1]],shift = 1:200)
  }
}  
s3 <- do.call(c,mclapply(cc3,noise_est,mc.cores = mc))




pdf(file = file.path(figs_dir,"rep1_npos_bet_50_100.pdf"),width = 6 ,height = 9)
dev.off()

## rep1, npos between 25 and 50
set.seed(123321)
keys <- stats[[1]][ between(npos,25,50),sample(match,500)]
cc4 <- list()

for(k in keys){
  if( stats[[1]][k, (label)] != "both"){
    cc4[[k]] <- NA
  }else{
    cc4[[k]] <- loc_scc_wrap(k,replicates[[1]],shift = 1:200)
  }
}

s4 <- do.call(c,mclapply(cc4,noise_est,mc.cores = mc))

          
pdf(file = file.path(figs_dir,"rep1_npos_bet_25_50.pdf"),width = 6 ,height = 9)
u <- mcmapply(plot_cc,cc,keys,MoreArgs = list(stats[[1]]),mc.cores = mc,SIMPLIFY = FALSE)
Z <- lapply(u,print)
dev.off()


## rep1, npos between 25 and 50
set.seed(123321)
keys <- stats[[1]][ between(npos,10,25),sample(match,500)]
cc5 <- list()
for(k in keys){
  if( stats[[1]][k, (label)] != "both"  | k == "chr18_200606"){
    cc5[[k]] <- NA
  }else{
    cc5[[k]] <- loc_scc_wrap(k,replicates[[1]],shift = 1:200)
  }
}

s5 <- do.call(c,mclapply(cc5,noise_est,mc.cores = mc))


set.seed(123321)
keys <- stats[[1]][ between(npos,5,10),sample(match,500)]
cc6 <- list()

for(k in keys){
  if( stats[[1]][k, (label)] != "both"){
    cc6[[k]] <- NA
  }else{
    cc6[[k]] <- loc_scc_wrap(k,replicates[[1]],shift = 1:200)
  }
}

s6 <- do.call(c,mclapply(cc6,noise_est,mc.cores = mc))



ss <-rbind(data.table(s = s1 , g = 1),data.table(s = s2,g = 2),data.table(s = s3,g = 3),
           data.table(s  = s4, g= 4),data.table(s = s5,g = 5),data.table(s = s6, g = 6))

         

ggplot(ss , aes(as.factor(g), s))+geom_boxplot()
dev.off()

