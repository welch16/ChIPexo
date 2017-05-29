#!/bin/sh

seqdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo/sequences
outdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo/fimo_peaks
motifdr=/p/keles/MEME/motif_databases

## 

db=JASPAR_CORE_2009.meme
motif=MA0148.1

motifdb=$motifdr/$db

fimo -oc $outdir/FoxA1_jaspar_Rep1_full --motif $motif $motifdb $seqdir/carroll_FoxA1_mouseliver_Rep1_exo_peak_sequences2.fna
fimo -oc $outdir/FoxA1_jaspar_Rep2_full --motif $motif $motifdb $seqdir/carroll_FoxA1_mouseliver_Rep2_exo_peak_sequences2.fna
fimo -oc $outdir/FoxA1_jaspar_Rep3_full --motif $motif $motifdb $seqdir/carroll_FoxA1_mouseliver_Rep3_exo_peak_sequences2.fna


# fimo -oc $outdir/ERR336935_FOXA1_$db_full --motif $motif $motifdb $seqdir/ERR336935_exo_peak_sequences2.fna
# fimo -oc $outdir/ERR336942_FOXA1_$db_full --motif $motif $motifdb $seqdir/ERR336942_exo_peak_sequences2.fna
# fimo -oc $outdir/ERR336956_FOXA1_$db_full --motif $motif $motifdb $seqdir/ERR336956_exo_peak_sequences2.fna

