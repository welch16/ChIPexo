
## PCR bottleneck report

### PCR bottleneck		



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
PBC <- function(file,cap = Inf)
{
  require(GenomicAlignments)
  require(spp)
  message(file)
  rr = readGAlignmentsFromBam(file,param = NULL)
  ss1 = subset(rr,subset =strand(rr)=="+")
  ss2 = subset(rr,subset = strand(rr)=="-")
  tab1 = table(start(ss1))
  tab2 = table(end(ss2))
  tab = table(c(start(ss1),end(ss2)))
  PBC1 = round(sum(tab1 == 1)/sum(tab1 >= 1 & tab1 <= cap),4)
  PBC2 = round(sum(tab2 == 1)/sum(tab2 >= 1 & tab2 <= cap),4)
  chip.data = read.bam.tags(file)
  table.chip.data = lapply(chip.data$tags,table)
  nUniq = sum(sapply(table.chip.data,function(x) sum(x==1)))
  nTotal =  sum(sapply(table.chip.data,length))
  PBC = round(nUniq/nTotal ,4)
  message("PBC +:",PBC1)
  message("PBC -:",PBC2)
  message("PBC:",PBC)
  return(c(PBC_plus=PBC1,PBC_minus=PBC2,PBC=PBC))  
}
```


### Figures

#### Figure 1

![plot of chunk overall_boxplot](figure/overall_boxplot.png) 

We can see:
- The overall level of the ChIP-Exo's PBC coefficients is lower than the level of the other two ChIP-Seq labels
- The level of ChIP-Seq PET is higher that the level of both SET generated y ChIP-Seq and ChIP-Exo

#### Figure 2

![plot of chunk boxplot_strand](figure/boxplot_strand.png) 

- For this plot, the PBC coefficents was calculated using only the reads with certain strand. In particular, we can see that for all 3 types, the boxplots for both strands are very similar.
- The observations of figure 1 hold for both strands.
- We can see that for all data types: The PBC level for each individual strand is very close to the one calculated using both strands.
- After inspecting the calculated values in the tables below, we can see that for both ChIP-Seq tags, the level of the PBC calculated with both strands is in between of the PBC's calculated using only one strand. However, for most of the ChIP-Exo samples, we can see that the PBC calculated with both strand is slightly higher than the PBC calculated for both individual strands.

### Tables

#### ChIP - Seq - SET

|                                        | PBC_plus| PBC_minus|    PBC|
|:---------------------------------------|--------:|---------:|------:|
|run101_Input_posO2_042814_qc.sorted.bam |   0.9773|    0.9763| 0.9768|
|run62_Input_negO2_042814_qc.sorted.bam  |   0.6431|    0.6434| 0.6433|
|run80_sig70_negO2_042814_qc.sorted.bam  |   0.5732|    0.5721| 0.5726|
|run80_sig70_posO2_042814_qc.sorted.bam  |   0.5611|    0.5600| 0.5606|

#### ChIP - Seq - PET

|                              | PBC_plus| PBC_minus|    PBC|
|:-----------------------------|--------:|---------:|------:|
|edsn1369_042814_qc.sorted.bam |   0.8194|    0.8192| 0.8193|
|edsn1396_042814_qc.sorted.bam |   0.7815|    0.7808| 0.7812|
|edsn1397_042814_qc.sorted.bam |   0.7369|    0.7359| 0.7364|
|edsn1398_042814_qc.sorted.bam |   0.6797|    0.6798| 0.6798|
|edsn1399_042814_qc.sorted.bam |   0.5086|    0.5089| 0.5087|
|edsn1400_042814_qc.sorted.bam |   0.7328|    0.7317| 0.7323|
|edsn1401_042814_qc.sorted.bam |   0.7367|    0.7361| 0.7364|
|edsn1402_042814_qc.sorted.bam |   0.7167|    0.7180| 0.7174|
|edsn1403_042814_qc.sorted.bam |   0.5490|    0.5493| 0.5492|
|edsn1416_042814_qc.sorted.bam |   0.8084|    0.8086| 0.8085|
|edsn788_042814_qc.sorted.bam  |   0.6599|    0.6597| 0.6598|
|edsn790_042814_qc.sorted.bam  |   0.4153|    0.4152| 0.4153|

#### ChIP - Exo

|                              | PBC_plus| PBC_minus|    PBC|
|:-----------------------------|--------:|---------:|------:|
|edsn1310_042814_qc.sorted.bam |   0.2849|    0.2805| 0.2904|
|edsn1311_042814_qc.sorted.bam |   0.2800|    0.2846| 0.2886|
|edsn1312_042814_qc.sorted.bam |   0.2545|    0.2482| 0.2561|
|edsn1313_042814_qc.sorted.bam |   0.3393|    0.3376| 0.3427|
|edsn1314_042814_qc.sorted.bam |   0.2691|    0.2706| 0.2749|
|edsn1315_042814_qc.sorted.bam |   0.2332|    0.2355| 0.2480|
|edsn1316_042814_qc.sorted.bam |   0.2587|    0.2529| 0.2635|
|edsn1317_042814_qc.sorted.bam |   0.2656|    0.2657| 0.2773|
|edsn1318_042814_qc.sorted.bam |   0.3036|    0.3011| 0.3048|
|edsn1319_042814_qc.sorted.bam |   0.4975|    0.4976| 0.4992|
|edsn1320_042814_qc.sorted.bam |   0.2151|    0.2157| 0.2169|
|edsn1321_042814_qc.sorted.bam |   0.2140|    0.2157| 0.2172|
|edsn930_042814_qc.sorted.bam  |   0.4964|    0.4985| 0.4969|
|edsn931_042814_qc.sorted.bam  |   0.2136|    0.2177| 0.2170|
|edsn932_042814_qc.sorted.bam  |   0.2950|    0.2939| 0.2995|
|edsn933_042814_qc.sorted.bam  |   0.2489|    0.2455| 0.2506|
|edsn934_042814_qc.sorted.bam  |   0.1749|    0.1783| 0.1799|
|edsn935_042814_qc.sorted.bam  |   0.2604|    0.2576| 0.2602|
|edsn936_042814_qc.sorted.bam  |   0.2323|    0.2307| 0.2380|
|edsn937_042814_qc.sorted.bam  |   0.2480|    0.2456| 0.2496|
|edsn938_042814_qc.sorted.bam  |   0.1633|    0.1621| 0.1679|
