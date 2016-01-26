#!/bin/sh

## performs chip exo analysis with MACE.

bamdir=/p/keles/ChIPexo/volume7/Landick/ChIPexo/rif_treatment
ecoli=inst/mace_analysis/ecoli.size
outdir=inst/mace_analysis


### rif-0min

rep1=$bamdir/edsn1311_Sig70.sort.bam
rep2=$bamdir/edsn1317_Sig70.sort.bam
outfile=Sig70_rif0

preprocessor.py -i $rep1,$rep2 -r $ecoli -o $outdir/$outfile

mace.py -s $ecoli -f $outdir/$outfile"_Forward.bw" -r $outdir/$outfile"_Reverse.bw" -o $outdir/$outfile

### rif-20min

rep1=$bamdir/edsn1314_Sig70.sort.bam
rep2=$bamdir/edsn1320_Sig70.sort.bam

outdir=inst/mace_analysis
outfile=Sig70_rif20

preprocessor.py -i $rep1,$rep2 -r $ecoli -o $outdir/$outfile

mace.py -s $ecoli -f $outdir/$outfile"_Forward.bw" -r $outdir/$outfile"_Reverse.bw" -o $outdir/$outfile
