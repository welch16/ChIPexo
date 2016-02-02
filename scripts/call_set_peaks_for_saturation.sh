#!/bin/sh

indir=/p/keles/ChIPexo/volume6/saturation_rif
fdr=.1
seed=2513

indir1=$indir/ChIPseq_SET/seed$seed
outdir1=$indir1/peaks

mkdir $outdir1
outdir1=$outdir1/FDR$fdr
mkdir $outdir1

rscripts/scripts/call_peaks.R $indir1 $outdir1 $fdr exo