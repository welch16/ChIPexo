#1/bin/sh

basedir=/p/keles/ChIPexo/volume4/tbp_analysis
seqexodir=$basedir/sequences/chipexo
seqnexusdir=$basedir/sequences/chipnexus

motifdr=/p/keles/MEME/motif_databases
db=JASPAR_CORE_2009.meme
motifdb=$motifdr/$db

fimoexodir=$basedir/fimo/chipexo
fimonexusdir=$basedir/fimo/chipnexus

## motif1

motif=MA0108.1

fimo -oc $fimoexodir/Rep1_TBP1 --motif $motif $motifdb $seqexodir/TBP_K562_Rep1.sequences.fna
fimo -oc $fimoexodir/Rep2_TBP1 --motif $motif $motifdb $seqexodir/TBP_K562_Rep2.sequences.fna
fimo -oc $fimoexodir/Rep3_TBP1 --motif $motif $motifdb $seqexodir/TBP_K562_Rep3.sequences.fna
fimo -oc $fimonexusdir/Rep1_TBP1 --motif $motif $motifdb $seqnexusdir/ChIPnexus_K562_TBP_rep1.sequences.fna
fimo -oc $fimonexusdir/Rep2_TBP1 --motif $motif $motifdb $seqnexusdir/ChIPnexus_K562_TBP_rep2.sequences.fna

## motif2

motif=MA0108.1

fimo -oc $fimoexodir/Rep1_TBP2 --motif $motif $motifdb $seqexodir/TBP_K562_Rep1.sequences.fna
fimo -oc $fimoexodir/Rep2_TBP2 --motif $motif $motifdb $seqexodir/TBP_K562_Rep2.sequences.fna
fimo -oc $fimoexodir/Rep3_TBP2 --motif $motif $motifdb $seqexodir/TBP_K562_Rep3.sequences.fna
fimo -oc $fimonexusdir/Rep1_TBP2 --motif $motif $motifdb $seqnexusdir/ChIPnexus_K562_TBP_rep1.sequences.fna
fimo -oc $fimonexusdir/Rep2_TBP2 --motif $motif $motifdb $seqnexusdir/ChIPnexus_K562_TBP_rep2.sequences.fna
