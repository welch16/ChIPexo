#!/bin/sh

# this script converts the new aligned reads to sam and generated the index bam files

basedir=/p/keles/ChIPexo/volume7/Landick
exodir1=$basedir/K12/ChIPexo/aerobic_vs_anaerobic
exodir2=$basedir/K12/ChIPexo/rif_treatment
petdir=$basedir/K12/ChIPseq_PET/rif_treatment

## old alignments ChIPexo

# scripts/format_bam.sh $exodir1/edsn931_Sig70.sam
# scripts/format_bam.sh $exodir1/edsn933_Sig70.sam
# scripts/format_bam.sh $exodir1/edsn935_Sig70.sam
# scripts/format_bam.sh $exodir1/edsn937_Sig70.sam

## new alignments ChIPexo

scripts/format_bam.sh $exodir2/edsn1311_Sig70.sam
scripts/format_bam.sh $exodir2/edsn1314_Sig70.sam
scripts/format_bam.sh $exodir2/edsn1317_Sig70.sam
scripts/format_bam.sh $exodir2/edsn1320_Sig70.sam

## new alignment ChIPseq_PET

scripts/format_bam.sh $petdir/edsn1369_Input.sam
scripts/format_bam.sh $petdir/edsn1396_Sig70.sam
scripts/format_bam.sh $petdir/edsn1398_Sig70.sam
scripts/format_bam.sh $petdir/edsn1400_Sig70.sam
scripts/format_bam.sh $petdir/edsn1402_Sig70.sam
