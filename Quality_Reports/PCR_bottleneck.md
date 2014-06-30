
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
And the following results where obtained:




```r
ls()
```

[1] "dr"     "EXO"    "f"      "files"  "folder" "PBC"    "PET"    "SET"   

