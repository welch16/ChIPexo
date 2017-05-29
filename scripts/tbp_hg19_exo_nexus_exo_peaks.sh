#!/bin/sh

basedir=/p/keles/ChIPexo/volume4/tbp_analysis
seqdir=$basedir/sequences

motifdr=/p/keles/MEME/motif_databases
db=JASPAR_CORE_2009.meme
motifdb=$motifdr/$db

seq1=$seqdir/TBP_K562_Rep1_exo_peak_sequences.fna
seq2=$seqdir/TBP_K562_Rep2_exo_peak_sequences.fna
seq3=$seqdir/TBP_K562_Rep3_exo_peak_sequences.fna
seq4=$seqdir/ChIPnexus_K562_TBP_rep1_exo_peak_sequences.fna
seq5=$seqdir/ChIPnexus_K562_TBP_rep2_exo_peak_sequences.fna

fimodir=$basedir/fimo_peaks

## motif1
motif=MA0108.1

fimo -oc $fimodir/venters_Rep1_jaspar --thresh 1e-2 --motif $motif $motifdb $seq1 &
fimo -oc $fimodir/venters_Rep2_jaspar --thresh 1e-2 --motif $motif $motifdb $seq2 &
fimo -oc $fimodir/venters_Rep3_jaspar --thresh 1e-2 --motif $motif $motifdb $seq3 &
fimo -oc $fimodir/nexus_Rep1_jaspar --thresh 1e-2 --motif $motif $motifdb $seq4 &
fimo -oc $fimodir/nexus_Rep2_jaspar --thresh 1e-2 --motif $motif $motifdb $seq5 &
