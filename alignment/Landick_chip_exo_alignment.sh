#!/bin/sh

index=$HOME/Desktop/Docs/Code/lib/bowtie2-2.2.6/indexes/ecoli/NC_000913
mc=24
maxAlign=1

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir1=$indir/run299
indir2=$indir/run423
outdir=$basedir/volume7
outdir1=$outdir/aerobic_vs_anaerobic
outdir2=$outdir/rif_treatment

# [welch@ramiz04] (113)$ bowtie2 --threads $mc -k $mismatch -x $index -U run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq.gz -S ex1320.sam

echo Alignining ChIP exo files:
echo Index file: $index
echo Threads used: $mc
echo Allowing reads to align to at most $maxAlign positions

## old files

# run299.456-166-1010-Sig-70-Lib_CAGATC_L007_R1.fastq.gz* 
# run299.456-166-1003b-Sig-70-Lib_GGCTAC_L007_R1.fastq.gz*
# run299.456-166-1003a-Sig-70-Lib_GATCAG_L007_R1.fastq.gz*
# run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq.gz* 

bowtie2 --threads $mc -k $maxAlign -x $index -U $indir1/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq.gz -S $outdir1/edsn931_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir1/run299.456-166-1010-Sig-70-Lib_CAGATC_L007_R1.fastq.gz -S $outdir1/edsn933_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir1/run299.456-166-1003a-Sig-70-Lib_GATCAG_L007_R1.fastq.gz -S $outdir1/edsn935_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir1/run299.456-166-1003b-Sig-70-Lib_GGCTAC_L007_R1.fastq.gz -S $outdir1/edsn937_Sig70.sam

## new files

# run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq.gz
# run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq.gz
# run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq.gz
# run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq.gz

bowtie2 --threads $mc -k $maxAlign -x $index -U $indir2/run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq.gz -S $outdir2/edsn1311_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir2/run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq.gz -S $outdir2/edsn1314_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir2/run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq.gz -S $outdir2/edsn1317_Sig70.sam
bowtie2 --threads $mc -k $maxAlign -x $index -U $indir2/run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq.gz -S $outdir2/edsn1320_Sig70.sam
