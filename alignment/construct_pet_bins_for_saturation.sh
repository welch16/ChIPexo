#!/bin/sh

indir=/p/keles/ChIPexo/volume6/saturation_rif
fragLen=150
binSize=150
seed=52314

indir2=$indir/ChIPseq_PET/seed$seed
outdir2=$indir2/bins
mkdir $outdir2
indir2=$indir2/samples

rscripts/construct_bins.R $indir2 $outdir2 $fragLen $binSize TRUE

