#!/bin/sh

basedir=/p/keles/ChIPexo/volume7/Landick
indir=$basedir/ChIPseq_PET/rif_treatment
outdir=$basedir/ChIPseq_SET/rif_treatment

rscripts/scripts/generate_SET.R $indir/edsn1369_Input.sort.bam $outdir/edsn1369_Input.bam
rscripts/scripts/generate_SET.R $indir/edsn1396_Sig70.sort.bam $outdir/edsn1396_Sig70.bam
rscripts/scripts/generate_SET.R $indir/edsn1398_Sig70.sort.bam $outdir/edsn1398_Sig70.bam
rscripts/scripts/generate_SET.R $indir/edsn1400_Sig70.sort.bam $outdir/edsn1400_Sig70.bam
rscripts/scripts/generate_SET.R $indir/edsn1402_Sig70.sort.bam $outdir/edsn1402_Sig70.bam

bamtools sort -in $outdir/edsn1369_Input.bam -out $outdir/edsn1369_Input.sort.bam
bamtools sort -in $outdir/edsn1396_Sig70.bam -out $outdir/edsn1396_Sig70.sort.bam
bamtools sort -in $outdir/edsn1398_Sig70.bam -out $outdir/edsn1398_Sig70.sort.bam
bamtools sort -in $outdir/edsn1400_Sig70.bam -out $outdir/edsn1400_Sig70.sort.bam
bamtools sort -in $outdir/edsn1402_Sig70.bam -out $outdir/edsn1402_Sig70.sort.bam

bamtools index -in $outdir/edsn1369_Input.sort.bam
bamtools index -in $outdir/edsn1396_Sig70.sort.bam
bamtools index -in $outdir/edsn1398_Sig70.sort.bam
bamtools index -in $outdir/edsn1400_Sig70.sort.bam
bamtools index -in $outdir/edsn1402_Sig70.sort.bam







