#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/saturation
fragLen=150
binSize=150
seed=$1

indir1=$indir/ChIPexo/seed$seed
outdir1=$indir1/bins
mkdir $outdir1
indir1=$indir1/samples

rscripts/scripts/construct_bins.R $indir1 $outdir1 $fragLen $binSize FALSE

