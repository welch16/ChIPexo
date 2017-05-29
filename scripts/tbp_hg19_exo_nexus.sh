#!/bin/sh

basedir=/p/keles/ChIPexo/volume4/tbp_analysis
seqdir=$basedir/sequences

motifdr=/p/keles/MEME/motif_databases
db=JASPAR_CORE_2009.meme
motifdb=$motifdr/$db

seq1=$seqdir/venters_TBP_K562_Rep1sequences.s
seq2=$seqdir/venters_TBP_K562_Rep2sequences.s
seq3=$seqdir/venters_TBP_K562_Rep3sequences.s
seq4=$seqdir/chipnexus_K562_TBP_Rep1sequences.s
seq5=$seqdir/chipnexus_K562_TBP_Rep2sequences.s

fimodir=$basedir/fimo

## motif1
motif=MA0108.1

fimo -oc $fimodir/venters_Rep1_jaspar --thresh 1e-2 --motif $motif $motifdb $seq1 &
fimo -oc $fimodir/venters_Rep2_jaspar --thresh 1e-2 --motif $motif $motifdb $seq2 &
fimo -oc $fimodir/venters_Rep3_jaspar --thresh 1e-2 --motif $motif $motifdb $seq3 &
fimo -oc $fimodir/nexus_Rep1_jaspar --thresh 1e-2 --motif $motif $motifdb $seq4 &
fimo -oc $fimodir/nexus_Rep2_jaspar --thresh 1e-2 --motif $motif $motifdb $seq5 &


