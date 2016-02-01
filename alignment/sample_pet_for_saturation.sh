#!/bin/sh

indir=/p/keles/ChIPexo/volume7/Landick/ChIPseq_PET/rif_treatment
outdir=/p/keles/ChIPexo/volume6/saturation_rif/ChIPseq_PET
minN=100000
maxN=900000
inc=100000
isPET=TRUE
seed=12452

outdir=$outdir/seed$seed
mkdir $outdir

outdir=$outdir/samples
mkdir $outdir

rscripts/scripts/sample_from_experiment.R $indir/edsn1396_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1398_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1400_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/scripts/sample_from_experiment.R $indir/edsn1402_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
