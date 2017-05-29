#!/bin/sh

basedir=/p/keles/ChIPexo/volume4/tbp_analysis
seqdir=$basedir/sequences

motifdr=/p/keles/MEME/motif_databases
db=JASPAR_CORE_2009.meme
motifdb=$motifdr/$db

seq1=$seqdir/venters_TBP_K562_Rep1_exo_peak_sequences2.fna
seq2=$seqdir/venters_TBP_K562_Rep2_exo_peak_sequences2.fna
seq3=$seqdir/venters_TBP_K562_Rep3_exo_peak_sequences2.fna
seq4=$seqdir/ChIPnexus_TBP_K562_Rep1_exo_peak_sequences2.fna
seq5=$seqdir/ChIPnexus_TBP_K562_Rep2_exo_peak_sequences2.fna

fimodir=$basedir/fimo_peaks

## motif1
motif=MA0108.1

fimo -oc $fimodir/venters_Rep1_jaspar2 --thresh 1e-2 --motif $motif $motifdb $seq1 &
fimo -oc $fimodir/venters_Rep2_jaspar2 --thresh 1e-2 --motif $motif $motifdb $seq2 &
fimo -oc $fimodir/venters_Rep3_jaspar2 --thresh 1e-2 --motif $motif $motifdb $seq3 &
fimo -oc $fimodir/nexus_Rep1_jaspar2 --thresh 1e-2 --motif $motif $motifdb $seq4 &
fimo -oc $fimodir/nexus_Rep2_jaspar2 --thresh 1e-2 --motif $motif $motifdb $seq5 &
