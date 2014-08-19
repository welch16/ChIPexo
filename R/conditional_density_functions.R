

#' @title A function to obtain the bin-counts qth quantile
#' @description This function divides the genome into bins and counts the overlaps between the bins and the reads in seq.set. Then return the probs-th quantiles of the counts
#' @param seq.set GRanges or GAlignments object with the reads obtained in a ChIP-Seq or ChIP-Exo experiment
#' @param binSize Numeric width of the bins
#' @param probs The probabilities for which the quantiles are going to be obtained
#' @return A numeric vector with the values of the quantiles
#' @export
seq.quantile <- function(seq.set,binSize,probs)
{  
  bins = create.bins(binSize,seqlengths(seq.set))
  counts = countOverlaps(bins,seq.set)
  return(quantile(counts,probs))  
}

#' @title A function to obtain the density of bin-counts given that the data is greater than it's prob-th quantile
#' @description This function uses \code{\link{density.reads.per.strand.ratio}}, where the reads are conditioned when the bin counts are greater than the prob-th quantile
#' @param binSize Numeric width of the bins
#' @param prob Numeric, the probability for which the quantile is obtained
#' @param type Character, the sequencing used to obtained the seqset
#' @param Rep Numeric, replicate number of data
#' @param seqset GRanges object with the aligned reads
#' @return A data.frame with the conditional density of the forward strand ratio and other characteristics
#' @export
quantile.df <- function(binSize,prob,type,Rep ,seqset)
{
  bins = create.bins(binSize,seqlengths(seqset))
  bins$counts = countOverlaps(bins,seqset)
  seqquantile = quantile(bins$counts,prob)
  dens = density.reads.per.strand.ratio(subset(bins,subset = counts > seqquantile),
    seqset)
  df = data.frame("Fwd.Strand.Ratio"=dens$x,density =dens$y,binSize = binSize,quantile = prob,Rep=Rep,type = type)
  return(df)
}

#' @title  This function repeats \code{\link{quantile.df}} to all values in binsizes and probs vectors
#' @param binsizes Numeric vector with all the desired binsizes
#' @param probs Numeric vector with all the desired probabilites
#' @param seqsets List of GRanges objects with a collection of reads
#' @param type Character, the sequencing used to obtained the seqset
#' @return A data frame
#' @seealso \code{\link{quantile.df}}
#' @export
quantile.multi.df <- function(binsizes,probs,seqsets,type,mc=8)
{
  
  all.binsizes = rep(binsizes,times = length(probs))
  all.prob = rep(probs,each =length(binsizes))

  seqlist = mclapply(1:length(seqsets),function(i,seqsets,all.binsizes,all.prob,type){
    df.list = mapply(quantile.df,all.binsizes,all.prob,MoreArgs = list(type,i,seqsets[[i]]),SIMPLIFY = FALSE)
    return(do.call(rbind,df.list))},seqsets,all.binsizes,all.prob,type,mc.cores = mc )

  df = do.call(rbind,seqlist)
  
  return(df)
}
