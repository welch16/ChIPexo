#!/bin/sh

indir=/p/keles/ChIPexo/volume7/Landick/ChIPseq_SET/rif_treatment
outdir=/p/keles/ChIPexo/volume6/saturation_rif/ChIPseq_SET
minN=100000
maxN=900000
inc=100000
isPET=FALSE
seed=12345

outdir=$outdir/seed$seed
mkdir $outdir

outdir=$outdir/samples
mkdir $outdir

rscripts/sample_from_experiment.R $indir/edsn1396_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment.R $indir/edsn1398_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment.R $indir/edsn1400_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
rscripts/sample_from_experiment.R $indir/edsn1402_Sig70.sort.bam $outdir $minN $maxN $inc $isPET $seed
