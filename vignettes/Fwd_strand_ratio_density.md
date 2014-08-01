
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

Initially several bin sizes were considered, but some of them undersmoothed the densities and others oversmoothed it. Therefore the densities calculated with bin size 200, 500 and 750 are shown. I

** Bin size = 200 **

```r
tab1 = subset(sample.info,eval(parse(text = st[[1]])))
edsn = as.character(tab1$edsn)
exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
exo.sets = lapply(exo.sets,function(y,exo)exo[[y]],exo)
pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
pet.sets = lapply(pet.sets,function(y,pet)pet[[y]],pet)
p1 = plot.density(200,exo.sets,pet.sets,genomeLength = seqlengths(exo.sets[[1]]))
print(p1)
```

![plot of chunk fig1](figure/fig1.png) 

** Bin size = 500 **

```r
tab1 = subset(sample.info,eval(parse(text = st[[1]])))
edsn = as.character(tab1$edsn)
exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
exo.sets = lapply(exo.sets,function(y,exo)exo[[y]],exo)
pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
pet.sets = lapply(pet.sets,function(y,pet)pet[[y]],pet)
p1 = plot.density(500,exo.sets,pet.sets,genomeLength = seqlengths(exo.sets[[1]]))
print(p1)
```

![plot of chunk fig2](figure/fig2.png) 

** Bin size = 750 **

```r
tab1 = subset(sample.info,eval(parse(text = st[[1]])))
edsn = as.character(tab1$edsn)
exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
exo.sets = lapply(exo.sets,function(y,exo)exo[[y]],exo)
pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
pet.sets = lapply(pet.sets,function(y,pet)pet[[y]],pet)
p1 = plot.density(750,exo.sets,pet.sets,genomeLength = seqlengths(exo.sets[[1]]))
print(p1)
```

![plot of chunk fig3](figure/fig3.png) 




#### Tables


|edsn |cult |ip            |phase       |growth  |rif    |rep |seq |
|:----|:----|:-------------|:-----------|:-------|:------|:---|:---|
|1311 |1197 |Sig70         |Exponential |Aerobic |0 min  |1   |Exo |
|1312 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min  |1   |Exo |
|1314 |1197 |Sig70         |Exponential |Aerobic |20 min |1   |Exo |
|1315 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |Exo |
|1317 |1202 |Sig70         |Exponential |Aerobic |0 min  |2   |Exo |
|1318 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min  |2   |Exo |
|1320 |1202 |Sig70         |Exponential |Aerobic |20 min |2   |Exo |
|1321 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |Exo |



|edsn |cult |ip            |phase       |growth  |rif    |rep |seq |
|:----|:----|:-------------|:-----------|:-------|:------|:---|:---|
|1396 |1197 |Sig70         |Exponential |Aerobic |0 min  |1   |PET |
|1397 |1197 |BetaPrimeFlag |Exponential |Aerobic |0 min  |1   |PET |
|1398 |1197 |Sig70         |Exponential |Aerobic |20 min |1   |PET |
|1399 |1197 |BetaPrimeFlag |Exponential |Aerobic |20 min |1   |PET |
|1400 |1202 |Sig70         |Exponential |Aerobic |0 min  |2   |PET |
|1401 |1202 |BetaPrimeFlag |Exponential |Aerobic |0 min  |2   |PET |
|1402 |1202 |Sig70         |Exponential |Aerobic |20 min |2   |PET |
|1403 |1202 |BetaPrimeFlag |Exponential |Aerobic |20 min |2   |PET |






