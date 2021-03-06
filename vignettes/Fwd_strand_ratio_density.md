
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{Forward Strand ratio density}
%\VignetteDepends{ggplot2, ChIP-Exo, GenomicAlignments}
-->

### Forward strand ratio density

For this analysis, we partitioned the E. Coli genome into bins. For each bin, the number of forward and backward reads that overlap it were counted. Finally, the forward strand proportion was calculated as:





#### Code used

Lets notice, a normal kernel is used and the bandwidth is estimated by cross - validation.



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

For this experiment we are comparing the samples when :
- The IP is either Sig70 or BetaPrimeFlag
- The samples were submitted to 0 min or 20 min of rif treatment
- The growth was aerobic
- The phase was exponential



#### Estimated densities

Initially several bin sizes were considered, but some of them undersmoothed the densities and others oversmoothed it. Therefore the densities calculated with bin size 200, 500 and 750 are shown. In particular we can see the following behaviour in allfigures:

- The densities for ChIP-Exo data sets seems to be almost uniform, specially with larger bin sizes
- The densities for ChIP-Exo and ChIP-Seq PET seems to be symmetrical, so it is not likely to be some sort of strand bias
- For ChIP-Exo, there are bins that overlap with a large amount of one strand reads but not the other. This may suggest some sort of enzyme digestion bias
- However, this is not related to a particular strand. We can asses that because for ChIP-Exo densities both tails are heavy and they seem to be of a very similar height





|edsn |cult |ip    |phase       |growth  |rif   |rep |seq |
|:----|:----|:-----|:-----------|:-------|:-----|:---|:---|
|1311 |1197 |Sig70 |Exponential |Aerobic |0 min |1   |Exo |
|1317 |1202 |Sig70 |Exponential |Aerobic |0 min |2   |Exo |
|1396 |1197 |Sig70 |Exponential |Aerobic |0 min |1   |PET |
|1400 |1202 |Sig70 |Exponential |Aerobic |0 min |2   |PET |

![plot of chunk fig1density](figure/fig1density.png) 




|edsn |cult |ip    |phase       |growth  |rif    |rep |seq |
|:----|:----|:-----|:-----------|:-------|:------|:---|:---|
|1314 |1197 |Sig70 |Exponential |Aerobic |20 min |1   |Exo |
|1320 |1202 |Sig70 |Exponential |Aerobic |20 min |2   |Exo |
|1398 |1197 |Sig70 |Exponential |Aerobic |20 min |1   |PET |
|1402 |1202 |Sig70 |Exponential |Aerobic |20 min |2   |PET |

![plot of chunk fig2density](figure/fig2density.png) 




|edsn |cult |ip            |phase       |growth  |rif   |rep |seq |
|:----|:----|:-------------|:-----------|:-------|:-----|:---|:---|
|1312 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min |1   |Exo |
|1318 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min |2   |Exo |
|1397 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min |1   |PET |
|1401 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min |2   |PET |

![plot of chunk fig3density](figure/fig3density.png) 




|edsn |cult |ip            |phase       |growth  |rif    |rep |seq |
|:----|:----|:-------------|:-----------|:-------|:------|:---|:---|
|1315 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |Exo |
|1321 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |Exo |
|1399 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |PET |
|1403 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |PET |

![plot of chunk fig4density](figure/fig4density.png) 





