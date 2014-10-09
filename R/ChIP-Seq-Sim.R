library(GenomicAlignments)

# The idea of this is to simulate a Chip - exo peak, to check when the strand cross correlation works as a method to find fragment length

peaksim <- function(nreads,mu,p,delta,sigma2,readLength)
{  
  strand = runif(nreads) >  p  # TRUE means Fwd, FALSE means Bwd
  fwd_pos = sum(strand)
  bwd_pos = nreads - fwd_pos
  fwd_start = round(rnorm(fwd_pos,mu - delta, sd = sqrt(sigma2)),0)
  bwd_start = round(rnorm(bwd_pos,mu + delta, sd = sqrt(sigma2)),0)
  peakReads = c(
    GRanges(seqnames = "sim",ranges = IRanges(start = fwd_start,width = readLength),strand = "+"),
    GRanges(seqnames = "sim",ranges = IRanges(end = bwd_start,width = readLength),strand = "-"))
  return(peakReads)
}

plotPeak <- function(peak,m)
{
  fwd_reads = suppressWarnings(subset(peak,subset = strand(peak)=="+"))
  bwd_reads = suppressWarnings(subset(peak,subset = strand(peak)=="-"))
  fwd_cov = coverage(fwd_reads)[[1]]
  bwd_cov = coverage(bwd_reads)[[1]]  
  x = 1:m
  fwd_xp = runLength(fwd_cov)[1:(nrun(fwd_cov))]
  fwd_yp = c(runValue(fwd_cov) ,0)
  fwd_y = stepfun(cumsum(fwd_xp),fwd_yp)(x)
  bwd_xp = runLength(bwd_cov)[1:(nrun(bwd_cov))]
  bwd_yp = c(runValue(bwd_cov),0)
  bwd_y = stepfun(cumsum(bwd_xp),bwd_yp)(x)
  mm = max(c(fwd_y,bwd_y))
  yl = 1.2 * c(0,mm)
  plot(NA,NA,ylim = yl,xlim = c(1,m),xlab = "genome cood.",ylab = "counts")
  lines(x,fwd_y,col = "red")
  lines(x,bwd_y,col = "blue")  
}  

# One iteration:

m = 1000  # Length of the genome
p = 0.5
delta = 50
sigma2 = 1000
nreads = 2000

mu = 500 #sample(delta:(m - delta) , 1   )

rl = 51
  

peak = peaksim(nreads = nreads,mu = mu , p = p,delta = delta,
               sigma2 = sigma2,readLength = rl)
layout(matrix(1:3,nrow = 3))
plotPeak(peak,m)
abline(v = mu,lty = 2)
cc =shiftApply(seq(1,300),
               coverage(subset(peak,subset = strand(peak) == "+"),width = m)[[1]],
               coverage(subset(peak,subset = strand(peak) == "-"),width = m)[[1]],
               cor,verbose = FALSE)
mm = which.max(cc)
plot(1:300,cc,type = "l",xlab="shift",ylab = "strand cc",main = mm)
plotPeak(resize(peak,mm),m)
abline(v = mu,lty = 2)


