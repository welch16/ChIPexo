#!/bin/sh

# this script converts the new aligned reads to sam and generated the index bam files

basedir=/p/keles/ChIPexo/volume7
dir1=$basedir/aerobic_vs_anaerobic
dir2=$basedir/rif_treatment

samtools view -b $dir1/edsn931_Sig70.sam > $dir1/edsn931_Sig70.bam
samtools view -b $dir1/edsn933_Sig70.sam > $dir1/edsn933_Sig70.bam
samtools view -b $dir1/edsn935_Sig70.sam > $dir1/edsn935_Sig70.bam
samtools view -b $dir1/edsn937_Sig70.sam > $dir1/edsn937_Sig70.bam

bamtools index -in $dir1/edsn931_Sig70.bam
bamtools index -in $dir1/edsn933_Sig70.bam
bamtools index -in $dir1/edsn935_Sig70.bam
bamtools index -in $dir1/edsn937_Sig70.bam

samtools view -b $dir2/edsn1311_Sig70.sam > $dir2/edsn1311_Sig70.bam
samtools view -b $dir2/edsn1314_Sig70.sam > $dir2/edsn1314_Sig70.bam
samtools view -b $dir2/edsn1317_Sig70.sam > $dir2/edsn1317_Sig70.bam
samtools view -b $dir2/edsn1320_Sig70.sam > $dir2/edsn1320_Sig70.bam
samtools view -b $dir2/edsn1369_Input.sam > $dir2/edsn1369_Input.bam
samtools view -b $dir2/edsn1396_Sig70.sam > $dir2/edsn1396_Sig70.bam
samtools view -b $dir2/edsn1398_Sig70.sam > $dir2/edsn1398_Sig70.bam
samtools view -b $dir2/edsn1400_Sig70.sam > $dir2/edsn1400_Sig70.bam
samtools view -b $dir2/edsn1402_Sig70.sam > $dir2/edsn1402_Sig70.bam

bamtools index -in $dir2/edsn1311_Sig70.bam
bamtools index -in $dir2/edsn1314_Sig70.bam
bamtools index -in $dir2/edsn1317_Sig70.bam
bamtools index -in $dir2/edsn1320_Sig70.bam
bamtools index -in $dir2/edsn1369_Input.bam
bamtools index -in $dir2/edsn1396_Sig70.bam
bamtools index -in $dir2/edsn1398_Sig70.bam
bamtools index -in $dir2/edsn1400_Sig70.bam
bamtools index -in $dir2/edsn1402_Sig70.bam
