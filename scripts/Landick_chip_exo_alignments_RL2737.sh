#!/bin/sh

index=/p/keles/ChIPexo/volume7/Landick/index/RL2737

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15/run423
outdir=$basedir/volume7/Landick/RL2737/rif_treatment/ChIPexo

mc=24

bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir/run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq $outdir/edsn1311_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir/run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq $outdir/edsn1314_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir/run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq $outdir/edsn1317_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir/run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq $outdir/edsn1320_Sig70.sam

