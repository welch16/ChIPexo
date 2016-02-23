#!/bin/sh

basedir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo

indir1=$basedir/rif_treatment
indir2=$basedir/aerobic_vs_anaerobic
outdir=/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo/bedfiles

bedtools bamtobed -i $indir1/edsn1311_Sig70.sort.bam > $outdir/edsn1311_Sig70.sort.bed
bedtools bamtobed -i $indir1/edsn1314_Sig70.sort.bam > $outdir/edsn1314_Sig70.sort.bed
bedtools bamtobed -i $indir1/edsn1317_Sig70.sort.bam > $outdir/edsn1317_Sig70.sort.bed
bedtools bamtobed -i $indir1/edsn1320_Sig70.sort.bam > $outdir/edsn1320_Sig70.sort.bed
bedtools bamtobed -i $indir2/edsn931_Sig70.sort.bam > $outdir/edsn931_Sig70.sort.bed
bedtools bamtobed -i $indir2/edsn933_Sig70.sort.bam > $outdir/edsn933_Sig70.sort.bed
