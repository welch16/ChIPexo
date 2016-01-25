#!/bin/sh

indir=/p/keles/ChIPexo/volume6/saturation_rif
fdr=.1
seed=12345

indir1=$indir/ChIPexo/seed$seed

readsdir=$indir1/samples

files=$(ls $readsdir/*bam) 

outdir1=$indir1/sites
mkdir $outdir1
outdir1=$outdir1/FDR$fdr
mkdir $outdir1

peaksdir=$indir1/peaks/FDR$fdr

for readsfile in $files
do
rscripts/call_binding_sites.R $readsfile $peaksdir $outdir1 FALSE 5
done

