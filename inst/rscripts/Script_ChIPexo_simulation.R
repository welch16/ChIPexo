
# We are going to follow dpeak's model with some small modifications to try to make ChIP exo experiments

rm(list = ls())
library(GenomicAlignments)
library(ggplot2)
library(data.table)

binding_sim <- function(nreads,mu,p,delta,sigma2,readlength)
{  
  strand = runif(nreads) >  p  # TRUE means Fwd, FALSE means Bwd
  fwd_pos = sum(strand)
  bwd_pos = nreads - fwd_pos
  if(fwd_pos > 0 ){
    fwd_start = round(rnorm(fwd_pos,mu - delta, sd = sqrt(sigma2)),0)
    fwd_reads = data.table(seqnames="sim",start=fwd_start,
      end=fwd_start+readlength-1,strand="+")
  }else{
    fwd_reads = data.table()
  }
  if(bwd_pos > 0){
    bwd_start = round(rnorm(bwd_pos,mu + delta, sd = sqrt(sigma2)),0)
    bwd_reads = data.table(seqnames="sim",start=bwd_start-readlength+1,
      end=bwd_start,strand ="-")
  }else{
    bwd_reads = data.table()
  }    
  peak_reads = rbind(fwd_reads,bwd_reads)
  return(peak_reads)
}


data.table2IRanges <- function(DT){
  ranges = IRanges(start = DT$start,
      end = DT$end)
  return(ranges)
}

data.table2GRanges <- function(DT){
  sn = as.character(DT$seqnames)
  gr = GRanges(seqnames = sn,
    ranges = data.table2IRanges(DT),
    strand = DT$strand)
  return(gr)
}


coverToVec <- function(lb,ub,cover)
{
  x = seq(lb,ub,by=1)
  xf = cumsum(runLength(cover))
  yf = c(runValue(cover),0)
  y = stepfun(xf,yf)(x)  
  return(y)
}





plotpeak <- function(peak,ext=100,mu=NULL)
{
  setkey(peak,strand)
  coords = c(peak$start,peak$end)
  lb = min(coords) - ext
  ub = max(coords) + ext
  fwd = peak["+"]
  bwd = peak["-"]
  vals = c()
  if(nrow(na.omit(fwd))>0){
    fwd_range = data.table2IRanges(fwd)
    fwd_cover = coverToVec(lb,ub,coverage(fwd_range))
    fwd_df = data.table(pos = lb:ub,counts=fwd_cover,strand="+")
    vals = c(vals,"+")
    m1 = max(fwd_df$counts)
  }else{
    fwd_df = data.table()
    m1 = 0
  }
  if(nrow(na.omit(bwd))>0){
    bwd_range = data.table2IRanges(bwd)
    bwd_cover = -coverToVec(lb,ub,coverage(bwd_range))
    bwd_df = data.table(pos = lb:ub,counts=bwd_cover,strand="-")
    vals = c(vals,"-")
    m2 = -min(bwd_df$counts)
  }else{
    bwd_df = data.table()
    m2 = 0
  }
  colors = c("blue","red")
  names(colors) = c("-","+")
  ylim = max(m1,m2)
  df = rbind(bwd_df,fwd_df)
  p = ggplot(df,aes(pos,counts,colour = strand))+
    geom_abline(intercept = 0,slope=0,linetype =2,size =.3)+
    geom_line()+theme_bw()+theme(legend.position ="none")+
      scale_colour_manual(values = colors[vals %in% names(colors)])+
      scale_y_continuous(limits = ylim*1.2*c(-1,1))
  if(!is.null(mu)){
    p = p + geom_vline(xintercept = mu,linetype = 3,size=.3)
  }
  return(p)
}




strand_statistics <- function(region_reads,param)
{
  # depth
  depth = nrow(region_reads)
  # Nr. positions
  ranges= data.table2IRanges(region_reads)
  if(param=="fwd"){
    nrPos = length(unique(start(ranges)))
  }else{
    nrPos = length(unique(end(ranges)))
  }
  # Summit position
  cover = coverage(ranges)
  if(nrun(cover) > 1){
    maxCover = max(cover)
    if(param=="fwd"){
      summitPos = head(which(cover == maxCover),n=1)
    }else{
      summitPos = tail(which(cover == maxCover),n=1)
    }
  }else{
    maxCover = NA
    summitPos = NA
  }
  return(c(depth=depth,nrPos =nrPos,maxCover=maxCover,summitPos=summitPos))
}

summary_statistics <- function(peak)
{
  setkey(peak,strand)
  fwd = peak["+"]
  bwd = peak["-"]
  fwd_stats = c(0,0,NA,NA)
  bwd_stats = c(0,0,NA,NA)
  if(nrow(na.omit(fwd))>0){
    fwd_stats = strand_statistics(fwd,"fwd")
  }
  if(nrow(na.omit(bwd))>0){
    bwd_stats = strand_statistics(bwd,"bwd")
  }
  depth = fwd_stats[1] + bwd_stats[1]
  nrPos = fwd_stats[2] + bwd_stats[2]
  prob = fwd_stats[1] / depth
  diff = bwd_stats[4] - fwd_stats[4]
  out = c(fwd_stats,bwd_stats,depth,nrPos,prob,diff)
  names(out) = c("f","f_nrPos","f_maxCover","f_summitPos",
         "r","r_nrPos","r_maxCover","r_summitPos",
         "depth","nrPos","prob","diff")         
  return(out)
}

p = 0.5
delta = 40
sigma2 = 1e4
nreads = 50
mu = 500
rl = 51
ext=10

peak = mclapply(1:1000,function(i)binding_sim(nreads = nreads,mu = mu , p = p,delta = delta,
  sigma2 = sigma2,readlength = rl),mc.cores =12)


plotpeaks <- function(...,ext=10,mu=NULL,mfrow = c(length(...),1))
{
  require(grid)
  mypos <- expand.grid(1:mfrow[1], 1:mfrow[2])
  mypos <- mypos[with(mypos, order(Var1)), ]
  pushViewport(viewport(layout = grid.layout(mfrow[1], mfrow[2])))
  
  plots = lapply(...,function(x)plotpeak(x,ext,mu))

  j <- 1
  for (i in 1:length(plots)){
    print(plots[[i]], vp=viewport(layout.pos.row=mypos[j,][1], layout.pos.col=mypos[j,][2]))
    j <- j+1
  }  
}

x11(width =12,height = 6)
plotpeaks(peak[sample(1000,4)],ext= 10,mu=500,mfrow = c(2,2))
                      



DT = data.table(do.call(rbind,mclapply(peak,summary_statistics,mc.cores=12)))
var(DT[,list(diff)][[1]])
mean(DT[,list(diff)][[1]])


