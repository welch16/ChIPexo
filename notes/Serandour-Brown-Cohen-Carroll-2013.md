
## Development of an Ilumina-based ChI-exonuclease method provides insight into FoxA1-DNA binding properties

The main resul is that ChIP-Exo outperforms ChIP-Seq in all the measures they choose. 

### Results

The authors adapted another ChIP-Seq method (different from Pugh's original metethod).

After applying MACS, the authors found:
* Almost all of the ChIP-Seq peaks are ChIP-Exo peaks.
* There is a significant amount of ChIP-Exo peaks that are not ChIP-Seq peaks
* The tails of the read depth densities of a ChIP-Seq experiment are heavier than the tails for a ChIP-Exo experiment
* Comparing the motif frequency around ChIP-Seq summits vs ChIP-Exo summits, the behaviour for the second one seemes to be more localized around the summit, sshowing a stronger signal around the summit (this is shown in figure 2 D-E) 

*Note - In this case they used the example of ER, FoxA1 and GATA3*

A challenge for ChIP-Seq experiment is to find a structure between functionally related transcription factors that bind to adjacent sequences and operate as acomplex.

There seems to be a sudden increase in read depth at a particular positions. The authors define this patters as 'mesas'. In a ChIP-Seq experiment this pattern is usually associated with an amplification artefact. In this example, there is a high frequency of motifs exactly 9 bp of the start position of the mesa, which suggest that is not an artefact.

* Also, its seems that the strand of the mesa is correlated with the orientation of the motif


### Note from authors

The main differences between this protocol and Pugh's protocol are the oligonucleotides sequences, different washing buffers, the use of magnetic bedsand the PCR mix.

The authors found no significant difference in peak width with increased exonuclease concentration.


[This paper could be downloaded here](http://genomebiology.com/content/pdf/gb-2013-14-12-r147.pdf)
