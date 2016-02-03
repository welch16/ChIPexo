#!/bin/sh

basedir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo

indir=$basedir/rif_treatment
outdir=$basedir/bedfiles

bedtools bamtobed -i $indir/edsn1311_Sig70.sort.bam > $outdir/edsn1311_Sig70.sort.bed
bedtools bamtobed -i $indir/edsn1314_Sig70.sort.bam > $outdir/edsn1314_Sig70.sort.bed
bedtools bamtobed -i $indir/edsn1317_Sig70.sort.bam > $outdir/edsn1317_Sig70.sort.bed
bedtools bamtobed -i $indir/edsn1320_Sig70.sort.bam > $outdir/edsn1320_Sig70.sort.bed

