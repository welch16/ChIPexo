#1/bin/sh

seqdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo/sequences
outdir=/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo
motifdr=/p/keles/MEME/motif_databases

## 

db=JASPAR_CORE_2009.meme
motif=MA0148.1

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA1_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna

db=JASPAR_CORE_2008.meme
motif=MA0047

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna

db=JASPAR_CORE_2009.meme
motif=MA0047.1

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna

db=JASPAR_CORE_2009.meme
motif=MA0047.2

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna

db=uniprobe_mouse.meme
motif=UP00073_1

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna

db=uniprobe_mouse.meme
motif=UP00073_2

motifdb=$motifdr/$db
fimo -oc $outdir/fimo/ERR336935_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336935.sequences.fna
fimo -oc $outdir/fimo/ERR336942_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336942.sequences.fna
fimo -oc $outdir/fimo/ERR336956_FOXA2_$db --motif $motif $motifdb $seqdir/ERR336956.sequences.fna
