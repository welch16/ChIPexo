#!/bin/sh

indir=/p/keles/ChIPexo/volume6/saturation_rif
fragLen=150
binSize=150
seed=12345

indir1=$indir/ChIPexo/seed$seed
outdir1=$indir1/bins
mkdir $outdir1
indir1=$indir1/samples

rscripts/construct_bins.R $indir1 $outdir1 $fragLen $binSize FALSE

indir2=$indir/ChIPseq_PET/seed$seed
outdir2=$indir2/bins
mkdir $outdir2
indir2=$indir2/samples

rscripts/construct_bins.R $indir2 $outdir2 $fragLen $binSize TRUE

indir3=$indir/ChIPseq_SET/seed$seed
outdir3=$indir3/bins
mkdir $outdir3
indir3=$indir3/samples

rscripts/construct_bins.R $indir3 $outdir3 $fragLen $binSize FALSE
