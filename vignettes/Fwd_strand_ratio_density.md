
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

For this experiment we are comparing the samples when:
- The IP is either Sig70 or BetaPrimeFlag
- The samples were submitted to 0 min or 20 min of rif treatment
- The growth was aerobic
- The phase was exponential




|   |edsn |cult |ip            |phase       |growth  |rif    |rep |seq |
|:--|:----|:----|:-------------|:-----------|:-------|:------|:---|:---|
|2  |1311 |1197 |Sig70         |Exponential |Aerobic |0 min  |1   |Exo |
|3  |1312 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min  |1   |Exo |
|5  |1314 |1197 |Sig70         |Exponential |Aerobic |20 min |1   |Exo |
|6  |1315 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |Exo |
|8  |1317 |1202 |Sig70         |Exponential |Aerobic |0 min  |2   |Exo |
|9  |1318 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min  |2   |Exo |
|11 |1320 |1202 |Sig70         |Exponential |Aerobic |20 min |2   |Exo |
|12 |1321 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |Exo |



|   |edsn |cult |ip            |phase       |growth  |rif    |rep |seq |
|:--|:----|:----|:-------------|:-----------|:-------|:------|:---|:---|
|23 |1396 |1197 |Sig70         |Exponential |Aerobic |0 min  |1   |PET |
|24 |1397 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min  |1   |PET |
|25 |1398 |1197 |Sig70         |Exponential |Aerobic |20 min |1   |PET |
|26 |1399 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |PET |
|27 |1400 |1202 |Sig70         |Exponential |Aerobic |0 min  |2   |PET |
|28 |1401 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min  |2   |PET |
|29 |1402 |1202 |Sig70         |Exponential |Aerobic |20 min |2   |PET |
|30 |1403 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |PET |








