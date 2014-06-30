
## PCR bottleneck report

### Index

[PCR bottleneck](#Def)


### PCR bottleneck [Def] ##



According to the **ENCODE** definitions of [quality metrics](http://encodeproject.org/ENCODE/qualityMetrics.html#definitions), one of the metrics used is the *PCR bottleneck coefficient* or **PBC**.

The definition is:

A measure of library complexity, i.e. how skewed the distribution of read counts per location is towards 1 read per location.

PBC = N1/Nd

(where N1= number of genomic locations to which **EXACTLY** one unique mapping read maps, and Nd = the number of genomic locations to which **AT LEAST** one unique mapping read maps, i.e. the number of non-redundant, unique mapping reads).

In particular we can see that always N1 <= Nd, so 0<= PBC <= 1. ENCODE recomends the following classification:

| PBC range | Bottleneck class |
| :---:     | :---: |
|0 -0.5  | Severe |
| 0.5-0.8 | Moderate|
|0.8-0.9 | Mild |
|0.9 - 1| Non-existant|


### Code

The following function was used to calculate the PBC:


```r
PBC <- function(file)
{
  require(GenomicAlignments)
  rr = readGAlignmentsFromBam(file,param = NULL)
  ss1 = subset(rr,subset =strand(rr)=="+")
  ss2 = subset(rr,subset = strand(rr)=="-")
  N11 = length(unique(start(ss1)))
  N12 = length(unique(end(ss2)))
  N1 = N11+N12
  Nd1 = length(ss1)
  Nd2 = length(ss2)
  Nd = Nd1 + Nd2
  PBC1 = round(N11 / Nd1,4)
  PBC2 = round(N12 / Nd2,4)
  PBC = round(N1 /Nd,4)
  return(c(PBC_plus=PBC1,PBC_minus=PBC2,PBC=PBC))  
}
```


### Figures

#### Figure 1

![plot of chunk overall_boxplot](figure/overall_boxplot.png) 

We can see:
- The overall level of the ChIP-Exo's PBC coefficients is lower than the level of the other two ChIP-Seq labels
- The range of the ChIP-Exo's PBC coefficients is also smaller than other two ranges. In particular, the interquartile range of the PET coefficients is smaller than the interquartile range of the SET coefficients.

#### Figure 2

![plot of chunk boxplot_strand](figure/boxplot_strand.png) 

- For this plot, the PBC coefficents was calculated using only the reads with certain strand. In particular, we can see that for all 3 types, the boxplots for both strands are very similar.
- The observations of figure 1 hold for both strands.



### Tables

#### ChIP - Seq - SET

|                                        | PBC_plus| PBC_minus|    PBC|
|:---------------------------------------|--------:|---------:|------:|
|run101_Input_posO2_042814_qc.sorted.bam |   0.9770|    0.9760| 0.9765|
|run62_Input_negO2_042814_qc.sorted.bam  |   0.6626|    0.6628| 0.6627|
|run80_sig70_negO2_042814_qc.sorted.bam  |   0.3125|    0.3131| 0.3128|
|run80_sig70_posO2_042814_qc.sorted.bam  |   0.3165|    0.3174| 0.3169|

#### ChIP - Seq - PET

|                              | PBC_plus| PBC_minus|    PBC|
|:-----------------------------|--------:|---------:|------:|
|edsn1369_042814_qc.sorted.bam |   0.8239|    0.8238| 0.8239|
|edsn1396_042814_qc.sorted.bam |   0.5972|    0.5977| 0.5975|
|edsn1397_042814_qc.sorted.bam |   0.5480|    0.5482| 0.5481|
|edsn1398_042814_qc.sorted.bam |   0.5222|    0.5227| 0.5225|
|edsn1399_042814_qc.sorted.bam |   0.2440|    0.2441| 0.2441|
|edsn1400_042814_qc.sorted.bam |   0.5205|    0.5205| 0.5205|
|edsn1401_042814_qc.sorted.bam |   0.5510|    0.5513| 0.5511|
|edsn1402_042814_qc.sorted.bam |   0.5719|    0.5727| 0.5723|
|edsn1403_042814_qc.sorted.bam |   0.3237|    0.3239| 0.3238|
|edsn1416_042814_qc.sorted.bam |   0.8139|    0.8141| 0.8140|
|edsn788_042814_qc.sorted.bam  |   0.3795|    0.3796| 0.3796|
|edsn790_042814_qc.sorted.bam  |   0.1844|    0.1843| 0.1843|

#### ChIP - Exo

|                              | PBC_plus| PBC_minus|    PBC|
|:-----------------------------|--------:|---------:|------:|
|edsn1310_042814_qc.sorted.bam |   0.0275|    0.0277| 0.0276|
|edsn1311_042814_qc.sorted.bam |   0.0902|    0.0919| 0.0911|
|edsn1312_042814_qc.sorted.bam |   0.0729|    0.0730| 0.0729|
|edsn1313_042814_qc.sorted.bam |   0.0176|    0.0185| 0.0180|
|edsn1314_042814_qc.sorted.bam |   0.1226|    0.1251| 0.1238|
|edsn1315_042814_qc.sorted.bam |   0.0347|    0.0345| 0.0346|
|edsn1316_042814_qc.sorted.bam |   0.0257|    0.0256| 0.0257|
|edsn1317_042814_qc.sorted.bam |   0.0457|    0.0475| 0.0466|
|edsn1318_042814_qc.sorted.bam |   0.1731|    0.1766| 0.1748|
|edsn1319_042814_qc.sorted.bam |   0.0243|    0.0274| 0.0258|
|edsn1320_042814_qc.sorted.bam |   0.0737|    0.0757| 0.0747|
|edsn1321_042814_qc.sorted.bam |   0.0421|    0.0435| 0.0428|
|edsn930_042814_qc.sorted.bam  |   0.0145|    0.0182| 0.0162|
|edsn931_042814_qc.sorted.bam  |   0.1019|    0.1024| 0.1022|
|edsn932_042814_qc.sorted.bam  |   0.0166|    0.0156| 0.0161|
|edsn933_042814_qc.sorted.bam  |   0.1111|    0.1076| 0.1093|
|edsn934_042814_qc.sorted.bam  |   0.0239|    0.0220| 0.0229|
|edsn935_042814_qc.sorted.bam  |   0.1298|    0.1292| 0.1295|
|edsn936_042814_qc.sorted.bam  |   0.0383|    0.0366| 0.0374|
|edsn937_042814_qc.sorted.bam  |   0.1124|    0.1119| 0.1121|
|edsn938_042814_qc.sorted.bam  |   0.0117|    0.0112| 0.0114|

