#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/saturation
fdr=.1
seed=$1

indir1=$indir/ChIPexo/seed$seed
outdir1=$indir1/peaks

mkdir $outdir1
outdir1=$outdir1/FDR$fdr
mkdir $outdir1

rscripts/scripts/call_peaks.R $indir1 $outdir1 $fdr exo
