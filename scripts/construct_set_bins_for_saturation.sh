#!/bin/sh

indir=/p/keles/ChIPexo/volume6/saturation_rif
fragLen=150
binSize=150
seed=2513

indir3=$indir/ChIPseq_SET/seed$seed
outdir3=$indir3/bins
mkdir $outdir3
indir3=$indir3/samples

rscripts/scripts/construct_bins.R $indir3 $outdir3 $fragLen $binSize FALSE
