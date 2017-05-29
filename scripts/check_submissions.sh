#!/bin/sh


basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir1=$indir/run299
indir2=$indir/run423

## first batch
md5sum $indir1/*Sig-70*fastq

# md5sum $indir1/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq
# md5sum $indir1/run299.456-166-1010-Sig-70-Lib_CAGATC_L007_R1.fastq 
# md5sum $indir1/run299.456-166-1003a-Sig-70-Lib_GATCAG_L007_R1.fastq
# md5sum $indir1/run299.456-166-1003b-Sig-70-Lib_GGCTAC_L007_R1.fastq

## second batch
md5sum $indir2/*ChIPexo*Sig*R1fastq

# md5sum $indir2/run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq
# md5sum $indir2/run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq
# md5sum $indir2/run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq
# md5sum $indir2/run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq

## wigfiles
wigdir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/wigfiles

md5sum $wigdir/*
