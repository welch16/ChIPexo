filterReads <- function(lb,ub,reads.pos)
  return(which(reads.pos >= lb & reads.pos <= ub))

filterSets <- function(exo.sets_cond,pet.sets_cond,distance)
{
  lengths = sapply(exo.sets_cond,length)
  ss = exo.sets_cond[[which.max(lengths)]]
  pos = start(ss)
  idx = which(diff(pos) > distance)

  lowerBounds = pmax(1,start(ss[idx])-distance)
  upperBounds = pmin(seqlengths(ss[idx]),start(ss[idx]) + distance)

  exo.idx = lapply(exo.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE))
  
  pet.idx = lapply(pet.sets_cond,function(y)
    mcmapply(FUN = filterReads,lowerBounds,upperBounds,MoreArgs = list(start(y)),SIMPLIFY = FALSE))

  return(list(exo = exo.idx,pet =pet.idx))  
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

single_set_plot <- function(lb,ub,reads,main="")
{
  reads_F = subset(reads,strand(reads)=="+")
  reads_R = subset(reads,strand(reads)=="-")
  window_f= coverToVec(lb,ub,coverage(reads_F)[[1]])
  window_r= coverToVec(lb,ub,coverage(reads_R)[[1]])
  x = seq(lb,ub,by=1)
  xlim = c(lb,ub)
  ylim = c(-1,1) * max(max(window_f),max(window_r))
  plot(x =1 ,y=0,xlim = xlim,ylim = ylim,type = "n",main =main ,xlab = "",ylab = "" )  
  abline(h=0,col =  "black",lty =2)
  lines(c(lb,x,ub),c(0,window_f,0),col = "red")
  lines(c(lb,x,ub),c(0,-window_r,0),col = "blue")  
}

all_plots <- function(lb,ub,exo.reads1,exo.reads2,pet.reads1,pet.reads2)
{
  par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
  single_set_plot(lb,ub,exo.reads1,"exo1")
  single_set_plot(lb,ub,exo.reads2,"exo2")
  single_set_plot(lb,ub,pet.reads1,"pet1")
  single_set_plot(lb,ub,pet.reads2,"pet2")
}
