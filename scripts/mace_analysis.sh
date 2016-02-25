#!/bin/sh

## performs chip exo analysis with MACE.

bamdir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo
bamdir1=$bamdir/rif_treatment
bamdir2=$bamdir/aerobic_vs_anaerobic
outdir=/p/keles/ChIPexo/volume6/K12/other_methods/mace
ecoli=$outdir/ecoli.size

## pre-process

preprocessor.py -i $bamdir1/edsn1311_Sig70.sort.bam -r $ecoli -o $outdir/edsn1311
preprocessor.py -i $bamdir1/edsn1314_Sig70.sort.bam -r $ecoli -o $outdir/edsn1314
preprocessor.py -i $bamdir1/edsn1317_Sig70.sort.bam -r $ecoli -o $outdir/edsn1317
preprocessor.py -i $bamdir1/edsn1320_Sig70.sort.bam -r $ecoli -o $outdir/edsn1320
preprocessor.py -i $bamdir2/edsn931_Sig70.sort.bam -r $ecoli -o $outdir/edsn931
preprocessor.py -i $bamdir2/edsn933_Sig70.sort.bam -r $ecoli -o $outdir/edsn933

## mace analysis

mace.py -s $ecoli -f $outdir/edsn1311_Forward.bw -r $outdir/edsn1311_Reverse.bw -o $outdir/edsn1311
mace.py -s $ecoli -f $outdir/edsn1314_Forward.bw -r $outdir/edsn1314_Reverse.bw -o $outdir/edsn1314
mace.py -s $ecoli -f $outdir/edsn1317_Forward.bw -r $outdir/edsn1317_Reverse.bw -o $outdir/edsn1317
mace.py -s $ecoli -f $outdir/edsn1320_Forward.bw -r $outdir/edsn1320_Reverse.bw -o $outdir/edsn1320
mace.py -s $ecoli -f $outdir/edsn931_Forward.bw -r $outdir/edsn931_Reverse.bw -o $outdir/edsn931
mace.py -s $ecoli -f $outdir/edsn933_Forward.bw -r $outdir/edsn933_Reverse.bw -o $outdir/edsn933
