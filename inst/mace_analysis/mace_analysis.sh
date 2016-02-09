#!/bin/sh

## performs chip exo analysis with MACE.

bamdir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment
ecoli=inst/mace_analysis/ecoli.size
outdir=inst/mace_analysis

## pre-process

preprocessor.py -i $bamdir/edsn1311_Sig70.sort.bam -r $ecoli -o $outdir/Sig70_rif0_rep1
preprocessor.py -i $bamdir/edsn1314_Sig70.sort.bam -r $ecoli -o $outdir/Sig70_rif20_rep1
preprocessor.py -i $bamdir/edsn1317_Sig70.sort.bam -r $ecoli -o $outdir/Sig70_rif0_rep2
preprocessor.py -i $bamdir/edsn1320_Sig70.sort.bam -r $ecoli -o $outdir/Sig70_rif20_rep2

## mace analysis

mace.py -s $ecoli -f $outdir/Sig70_rif0_rep1_Forward.bw -r $outdir/Sig70_rif0_rep1_Reverse.bw -o $outdir/Sig70_rif0_rep1
mace.py -s $ecoli -f $outdir/Sig70_rif20_rep1_Forward.bw -r $outdir/Sig70_rif20_rep1_Reverse.bw -o $outdir/Sig70_rif20_rep1
mace.py -s $ecoli -f $outdir/Sig70_rif0_rep2_Forward.bw -r $outdir/Sig70_rif0_rep2_Reverse.bw -o $outdir/Sig70_rif0_rep2
mace.py -s $ecoli -f $outdir/Sig70_rif20_rep2_Forward.bw -r $outdir/Sig70_rif20_rep2_Reverse.bw -o $outdir/Sig70_rif20_rep2
