

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

First we are showing a comparison of the conditional densities of the ChIP-seq pet samples vs the ChIP-exo samples such that their edsn is > 1300. Considering the case when the Ip is Sig70, the rif time is 20 min, the grown is aerobic and the phase is exponential we can see, only for the first replicate:








```
## Warning: Removed 758 rows containing missing values (geom_path).
## Warning: Removed 750 rows containing missing values (geom_path).
## Warning: Removed 739 rows containing missing values (geom_path).
## Warning: Removed 704 rows containing missing values (geom_path).
## Warning: Removed 248 rows containing missing values (geom_path).
```

![plot of chunk fig1_condDensity](figure/fig1_condDensity.png) 




