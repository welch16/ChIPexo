
## Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data

*Abstract* The idea is investigate the effects of blacklisting and duplicated reads on established metrics of ChIP-seq quality. There is also a first application of these metrics to ChIP-Exo data, and some recommendations.

### Blacklists:

The authors considered 3 blacklisted regions from where the first two can be downloaded from the UCSC genome browser
[here](http://hgwdev.cse.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability)
- DAC
- DER
- UHS [this can be downloaded from here](https://sites.google.com/site/anshulkundaje/projects/blacklists)

### Quality metrics and cross - correlation profiles:

- SSD metrics where calculated using [htSeqTools](http://www.bioconductor.org/packages/release/bioc/html/htSeqTools.html), this calculates dispersion measures dependend on the number of reads. Like the standard deviation or Gini's index. The idea is that by comparing this measures, inneficient inmunoprecipitation (due to antibody's quality) samples can be revealed. 
- RSC (Relative Strand Cross - Correlation) and NSC (Normalized Strand Cross - Correlation) where calculated using [ccQualityControl](https://code.google.com/p/phantompeakqualtools/), the definitions of both measures must be [here](http://www.nature.com/nbt/journal/v26/n12/full/nbt.1508.html) or [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3931556/)

### Results

**Def** Reads mapped to the same position are termed as *duplicates*

- Historically duplicated reads have been removed in order to remove artifacts from PCR amplification bias and sources of aberrant signal. This originated the NRF (Non - Redundant Fraction) measure, which suggest that less than 20% of the reads must be duplicates for 10 million reads sequenced.
- Blacklisting regions attempts to remove sources of artificial signal caused by biases of chromatin accesibility and ambiguous alignment.


#### Figures

- **Figure 2**: For A and B the reads are classified as all (red) (*need to check if it means all reads or the reads that aren't duplicated or multi-mapped*), duplicated (blue) and multimapped (orange):

    * A - In this figure, the authors showed that the reads for all the experiments analyzed in the blacklisted regions are either duplicated or multimapped.
  
    * B - This are boxplots of the same classification of reads, checking the percentage of total reads in DAC regions. We can see that if we compare ChIP-seq vs ChIP-Exo the quantity is smaller for ChIP-Exo in the 3 categories (all, duplicated and multimapped)
	
    * C - This are RPKM within blacklisted regions.

- **Figure 3**: In this plots the quantity of interest if the SSD

    * A - The SSD of ENCODE data is higher before being blacklisted

	* B - The order in variability of the SSD is Input  < TF  < Histone

- **Figure 4**: Illustration of extension of fragment length in both strands

- **Figure 5**: Example and explanation of the estimation of fragment length by the cross - correlation approach

- **Figure 6**: Scatter plots of the estimated fragment length before (the state is defined as the y axis) and after blacklisting (respect to x axis). For example, on the first panel the x-axis are the fragment length estimated after blacklisting with the DER blacklist and the y axis corresponds to the fragment lengths estimated using all the reads.



[The paper can be downloaded from here](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00075/full)
