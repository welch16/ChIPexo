#!/bin/sh

indir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment
outdir=/p/keles/ChIPexo/volume6/K12/saturation/ChIPexo
minN=100000
maxN=900000
inc=100000
isPET=FALSE
seed=$1

outdir=$outdir/seed$seed
mkdir $outdir

outdir=$outdir/samples
mkdir $outdir

rscripts/scripts/sample_from_experiment.R $indir/edsn1311_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1314_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1317_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1320_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
