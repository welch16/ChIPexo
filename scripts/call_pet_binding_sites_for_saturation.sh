#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/saturation
fdr=.1
seed=$1

indir1=$indir/ChIPseq_PET/seed$seed

readsdir=$indir1/samples

files=$(ls $readsdir/*bam) 

outdir1=$indir1/sites
mkdir $outdir1
outdir1=$outdir1/FDR$fdr
mkdir $outdir1

peaksdir=$indir1/peaks/FDR$fdr

for readsfile in $files
do
rscripts/scripts/call_binding_sites.R $readsfile $peaksdir $outdir1 TRUE 5
done
