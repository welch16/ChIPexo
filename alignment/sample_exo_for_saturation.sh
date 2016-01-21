#!/bin/sh

indir=/p/keles/ChIPexo/volume7/Landick/ChIPexo/rif_treatment
outdir=/p/keles/ChIPexo/volume6/saturation_rif/ChIPexo
minN=100000
maxN=900000
inc=50000
isPET=FALSE
seed=12345

outdir=$outdir/seed$seed
mkdir $outdir

rscripts/sample_from_experiment $infile/edsn1311_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment $infile/edsn1314_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment $infile/edsn1317_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment $infile/edsn1320_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
