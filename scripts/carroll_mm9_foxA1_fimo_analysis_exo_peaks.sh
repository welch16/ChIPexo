#!/bin/sh

seqdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo/sequences
outdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo/fimo_peaks
motifdr=/p/keles/MEME/motif_databases

## 

db=JASPAR_CORE_2009.meme
motif=MA0148.1

motifdb=$motifdr/$db
fimo -oc $outdir/ERR336935_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336935_exo_peak_sequences.fna
fimo -oc $outdir/ERR336942_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336942_exo_peak_sequences.fna
fimo -oc $outdir/ERR336956_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336956_exo_peak_sequences.fna

