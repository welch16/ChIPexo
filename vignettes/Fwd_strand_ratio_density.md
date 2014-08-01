
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{Forward Strand ratio density}
%\VignetteDepends{ggplot2, ChIP-Exo, GenomicAlignments}
-->

# Forward strand ratio density

For this analysis, we partitioned the E. Coli genome into bins. For each bin, the number of forward and backward reads that overlap it were counted. Finally, the forward strand proportion was calculated as:






```r
density.reads.per.strand.ratio
```

```
## function (bins, reads) 
## {
##     counts_F = countOverlaps(bins, subset(reads, subset = strand(reads) == 
##         "+"))
##     counts_R = countOverlaps(bins, subset(reads, subset = strand(reads) == 
##         "-"))
##     ratio = (counts_F + 1)/(counts_F + counts_R + 2)
##     return(density(ratio))
## }
```






