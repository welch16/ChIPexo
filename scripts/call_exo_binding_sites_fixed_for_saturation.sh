#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/saturation
edsn=$1
seed=$2
fdr=.1
indir1=$indir/ChIPexo/seed$seed

echo $indir1
echo $edsn

readsdir=$indir1/samples
files=$(ls $readsdir/edsn$edsn*bam) 

outdir1=$indir1/sites_fixed
mkdir $outdir1
outdir1=$outdir1/FDR$fdr
mkdir $outdir1

peaksdir=$indir1/peaks/FDR$fdr

peaksfile=$(ls $peaksdir/edsn$edsn*900K*)

for readsfile in $files
do
rscripts/scripts/call_binding_sites.R $readsfile $peaksdir $outdir1 FALSE 5
done
