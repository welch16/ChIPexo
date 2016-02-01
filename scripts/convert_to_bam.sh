#!/bin/sh

# this script converts the new aligned reads to sam and generated the index bam files

basedir=/p/keles/ChIPexo/volume7/Landick
exodir1=$basedir/ChIPexo/aerobic_vs_anaerobic
exodir2=$basedir/ChIPexo/rif_treatment
petdir=$basedir/ChIPseq_PET/rif_treatment

## old alignments ChIPexo

alignment/format_bam.sh $exodir1/edsn931_Sig70.sam
alignment/format_bam.sh $exodir1/edsn933_Sig70.sam
alignment/format_bam.sh $exodir1/edsn935_Sig70.sam
alignment/format_bam.sh $exodir1/edsn937_Sig70.sam

## new alignments ChIPexo

alignment/format_bam.sh $exodir2/edsn1311_Sig70.sam
alignment/format_bam.sh $exodir2/edsn1314_Sig70.sam
alignment/format_bam.sh $exodir2/edsn1317_Sig70.sam
alignment/format_bam.sh $exodir2/edsn1320_Sig70.sam

## new alignment ChIPseq_PET

alignment/format_bam.sh $petdir/edsn1369_Input.sam
alignment/format_bam.sh $petdir/edsn1396_Sig70.sam
alignment/format_bam.sh $petdir/edsn1398_Sig70.sam
alignment/format_bam.sh $petdir/edsn1400_Sig70.sam
alignment/format_bam.sh $petdir/edsn1402_Sig70.sam
