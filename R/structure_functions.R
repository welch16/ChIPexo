filterReads <- function(lb,ub,reads.pos)
  return(which(reads.pos >= lb & reads.pos <= ub))

filterSets <- function(exo.sets_cond,pet.sets_cond,distance,mc=8)
{
  lengths = sapply(exo.sets_cond,length)
  ss = exo.sets_cond[[which.max(lengths)]]
  pos = start(ss)
  idx = which(diff(pos) > distance)

  lowerBounds = pmax(1,start(ss[idx])-distance)
  upperBounds = pmin(seqlengths(ss[idx]),start(ss[idx]) + distance)

  exo.idx = lapply(exo.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE,mc.cores = mc))
  
  pet.idx = lapply(pet.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE,mc.cores = mc))

  return(list(exo = exo.idx,pet =pet.idx))  
}


filterSets2 <- function(exo.sets_cond,pet.sets_cond,set.sets_cond,distance,mc=8)
{
  lengths = sapply(exo.sets_cond,length)
  ss = exo.sets_cond[[which.max(lengths)]]
  pos = start(ss)
  idx = which(diff(pos) > distance)

  lowerBounds = pmax(1,start(ss[idx])-distance)
  upperBounds = pmin(seqlengths(ss[idx]),start(ss[idx]) + distance)

  exo.idx = lapply(exo.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE,mc.cores = mc))
  
  pet.idx = lapply(pet.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE,mc.cores = mc))

  set.idx = lapply(set.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE,mc.cores = mc))
  
  return(list(exo = exo.idx,pet =pet.idx,set=set.idx))  
}


getReads <- function(idx,reads,mc)mclapply(idx,function(x,reads)reads[x],reads,mc.cores =mc)
coverToVec <- function(lb,ub,cover)
{
  x = seq(lb,ub,by=1)
  if(nrun(cover)==1){
    y = rep(runValue(cover),length(x))
  }else{
    xf = cumsum(runLength(cover)[1:(nrun(cover)-1)])
    yf = runValue(cover)
    y = stepfun(xf,yf)(x)
  }
  return(y)
}

cover_by_strand <- function(reads,str){
  out <- suppressWarnings(coverage(subset(reads,strand(reads)==str))[[1]])
  return(out)
}
  
single_set_plot <- function(lb,ub,reads,ext,main="")
{
  if(is.null(reads)){
    x = seq(lb-ext,ub+ext,by=1)
    xlim = c(lb-ext,ub+ext)
    ylim = c(-1,1) 
    plot(x =1 ,y=0,xlim = xlim,ylim = ylim,type = "n",main =main ,xlab = "",ylab = "" )
    abline(h=0,col =  "black",lty =2)
    abline( v = c(lb,ub),lty=2)
  }else{
    cover_F = cover_by_strand(reads,"+")
    cover_R = cover_by_strand(reads,"-")
    window_f= coverToVec(lb-ext,ub+ext,cover_F)
    window_r= coverToVec(lb-ext,ub+ext,cover_R)
    x = seq(lb-ext,ub+ext,by=1)
    xlim = c(lb-ext,ub+ext)
    ylim = c(-1,1) * max(max(window_f),max(window_r))
    plot(x =1 ,y=0,xlim = xlim,ylim = ylim,type = "n",main =main ,
         xlab = "",ylab = "" )  
    abline(h=0,col =  "black",lty =2)
    lines(c(lb,x,ub),c(0,window_f,0),col = "red")
    lines(c(lb,x,ub),c(0,-window_r,0),col = "blue")
    abline(v = c(lb,ub),lty=2)
  }
}

all_plots <- function(lb,ub,exo.reads1,exo.reads2,pet.reads1,pet.reads2)
{
  par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
  single_set_plot(lb,ub,exo.reads1,"exo1")
  single_set_plot(lb,ub,exo.reads2,"exo2")
  single_set_plot(lb,ub,pet.reads1,"pet1")
  single_set_plot(lb,ub,pet.reads2,"pet2")
}


all_plots2 <- function(lb,ub,exo.reads1,exo.reads2,pet.reads1,pet.reads2,set.reads1,set.reads2,ext=0)
{
  par(mfcol=c(2,3), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
  single_set_plot(lb,ub,exo.reads1,ext,"exo1")
  single_set_plot(lb,ub,exo.reads2,ext,"exo2")
  single_set_plot(lb,ub,pet.reads1,ext,"pet1")
  single_set_plot(lb,ub,pet.reads2,ext,"pet2")
  single_set_plot(lb,ub,set.reads1,ext,"set1")
  single_set_plot(lb,ub,set.reads2,ext,"set2")
}
