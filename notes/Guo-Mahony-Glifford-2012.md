
## High Resolution Genome Wide Binding Event Finding and Motif Discovery Reveals Transcription Factor Spatial Binding Constraints

#### Note

For this paper, the authors developed [GEM](http://www.psrg.csail.mit.edu/gem/)

### Motivation

Genomic sequences facilitates both competitive and cooperative regulatory factor-factor interactions that implements transcriptional regulatory logic. Appropiately spaced motifs can facilitate cooperative binding, while overlapping motifs can implement competitive binding.

**Def** (Binding event location) Single base position at the center of a protein-DNA interaction.

### Results

**Def** (Spatial resolution) Average absolute value difference between predicted motif location and the nearest match to a proximal consensus motif

Claim: GEM has the best spatial resolution among the compared methods, the compared methods are older methods like MACS, cisGENOME, QuEST, etc. *Wonder how it would work vs mosaics + dpeak*

GEM accurately discovers DNA-binding motifs in ENCODE ChIP-Seq data

#### ChIP - Exo

GEM's model could be applied automatically to ChIP - Exo data (they applied it to the published CTCF data set), the distribution estimated by GEM agrees with the empirical distribution of the CTCF data sets.

### Methods:

An outline of GEM's method:

1. Predict protein-DNA binding event locations with a sparse prior
2. Discover the set of enriched k-mers at binding event locations
3. Cluster the set of enriched k-mers into k-mer equivalence classes
4. Generate a positional prior for event discovery with the most enriched k-mer equivalence class
5. Predict improved protein-DNA binding event locations with a k-mer based positional prioir
6. Repeat motif discovery (steps 2-3) from the Phase 5 improved event locations


[The paper could be downloaded from here](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002638)
