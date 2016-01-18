#!/bin/sh

sam=$1
bam=${1/.sam/.bam}

sortbam=${bam/.bam/.sort.bam}

samtools view -b $sam > $bam
bamtools sort -in $bam -out $sortbam
bamtools index -in $sortbam

