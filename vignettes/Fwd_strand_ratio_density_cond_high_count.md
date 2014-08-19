

<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{Forward Strand ratio density}
%\VignetteDepends{ggplot2, ChIP-Exo, GenomicAlignments}
-->

### Forward strand ratio density conditional to high count regions

For this analysis, we followed the same ideas as in the [forward strand ratio analysis](../Fwd_strand_ratio_density.md). However, we are conditioning the bins to consider only the ones that have counts higher than certain quantile.







```r
quantile.df
```

```
## function (binSize, prob, type, Rep, seqset) 
## {
##     bins = create.bins(binSize, seqlengths(seqset))
##     bins$counts = countOverlaps(bins, seqset)
##     seqquantile = quantile(bins$counts, prob)
##     dens = density.reads.per.strand.ratio(subset(bins, subset = counts > 
##         seqquantile), seqset)
##     df = data.frame(Fwd.Strand.Ratio = dens$x, density = dens$y, 
##         binSize = binSize, quantile = prob, Rep = Rep, type = type)
##     return(df)
## }
```



