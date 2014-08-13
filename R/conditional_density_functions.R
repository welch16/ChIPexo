

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

