#!/bin/sh

indir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic
outdir=/p/keles/ChIPexo/volume6/K12/saturation_QC/aerobic_vs_anaerobic
minN=10000
maxN=90000
inc=10000
isPET=FALSE
seed=$1

outdir=$outdir/seed$seed
mkdir $outdir

outdir=$outdir/samples
mkdir $outdir

rscripts/scripts/sample_from_experiment.R $indir/edsn931_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn933_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
# rscripts/scripts/sample_from_experiment.R $indir/edsn935_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
# rscripts/scripts/sample_from_experiment.R $indir/edsn037_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
