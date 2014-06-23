
## Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data

*Abstract* The idea is investigate the effects of blacklisting and duplicated reads on established metrics of ChIP-seq quality. There is also a first application of these metrics to ChIP-Exo data, and some recommendations.

### Blacklists:

The authors considered 3 blacklisted regions from where the first two can be downloaded from the UCSC genome browser
[here](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability)
- DAC
- DER
- UHS [this can be downloaded from here](https://sites.google.com/site/anshulkundaje/projects/blacklists)

### Quality metrics and cross - correlation profiles

- SSD metrics where calculated using [htSeqTools](http://www.bioconductor.org/packages/release/bioc/html/htSeqTools.html), this calculates dispersion measures dependend on the number of reads. Like the standard deviation or Gini's index. The idea is that by comparing this measures, inneficient inmunoprecipitation (due to antibody's quality) samples can be revealed. 
- RSC (Relative Strand Cross - Correlation) and NSC (Normalized Strand Cross - Correlation) where calculated using [ccQualityControl](https://code.google.com/p/phantompeakqualtools/), the definitions of both measures must be [here](http://www.nature.com/nbt/journal/v26/n12/full/nbt.1508.html) or [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3931556/)




[The paper can be downloaded from here](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00075/full)
