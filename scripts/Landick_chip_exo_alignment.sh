#!/bin/sh

index=/p/keles/ChIPexo/volume7/Landick/index/E.Coli_K-12
mc=24

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir1=$indir/run299
indir2=$indir/run423
outdir=$basedir/volume7/Landick/K12/ChIPexo
outdir1=$outdir/aerobic_vs_anaerobic
outdir2=$outdir/rif_treatment

# [welch@ramiz04] (113)$ bowtie2 --threads $mc -k $mismatch -x $index -U run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq.gz -S ex1320.sam

echo Alignining ChIP exo files:
echo Index file: $index
echo Threads used: $mc

## old files

# run299.456-166-1010-Sig-70-Lib_CAGATC_L007_R1.fastq.gz* 
# run299.456-166-1003b-Sig-70-Lib_GGCTAC_L007_R1.fastq.gz*
# run299.456-166-1003a-Sig-70-Lib_GATCAG_L007_R1.fastq.gz*
# run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq.gz* 
## -m 1 -l 55 -k 1 -5 3 -3 40 --best -S



bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq $outdir1/edsn931_Sig70.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1010-Sig-70-Lib_CAGATC_L007_R1.fastq $outdir1/edsn933_Sig70.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1003a-Sig-70-Lib_GATCAG_L007_R1.fastq $outdir1/edsn935_Sig70.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1003b-Sig-70-Lib_GGCTAC_L007_R1.fastq $outdir1/edsn937_Sig70.sam

# reads processed: 16846671
# reads with at least one reported alignment: 13961493 (82.87%)
# reads that failed to align: 2164204 (12.85%)
# reads with alignments suppressed due to -m: 720974 (4.28%)
# Reported 13961493 alignments to 1 output stream(s)
# reads processed: 17268284
# reads with at least one reported alignment: 14810838 (85.77%)
# reads that failed to align: 1788229 (10.36%)
# reads with alignments suppressed due to -m: 669217 (3.88%)
# Reported 14810838 alignments to 1 output stream(s)
# reads processed: 18868657
# reads with at least one reported alignment: 16108774 (85.37%)
# reads that failed to align: 1622603 (8.60%)
# reads with alignments suppressed due to -m: 1137280 (6.03%)
# Reported 16108774 alignments to 1 output stream(s)
# reads processed: 16423038
# reads with at least one reported alignment: 13636541 (83.03%)
# reads that failed to align: 1795625 (10.93%)
# reads with alignments suppressed due to -m: 990872 (6.03%)
# Reported 13636541 alignments to 1 output stream(s)


## new files

# run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq.gz
# run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq.gz
# run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq.gz
# run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq.gz

# -q -v 2 -a -m 1

bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq $outdir2/edsn1311_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1314-cult-1197-ChIPexo-20minRif-Sigma70_ACAGTG_L001_R1.fastq $outdir2/edsn1314_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1317-cult-1202-ChIPexo-0minRif-Sigma70_ACTTGA_L001_R1.fastq $outdir2/edsn1317_Sig70.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq $outdir2/edsn1320_Sig70.sam


# Index file: /p/keles/ChIPexo/volume7/Landick/index/E.Coli_K-12 
# Threads used: 24
# 1311
# reads processed: 1050680
# reads with at least one reported alignment: 902921 (85.94%)
# reads that failed to align: 89270 (8.50%)
# reads with alignments suppressed due to -m: 58489 (5.57%)
# Reported 902921 alignments to 1 output stream(s)
# 1314
# reads processed: 2115858
# reads with at least one reported alignment: 1852124 (87.54%)
# reads that failed to align: 188420 (8.91%)
# reads with alignments suppressed due to -m: 75314 (3.56%)
# Reported 1852124 alignments to 1 output stream(s)
# 1317
# reads processed: 2444898
# reads with at least one reported alignment: 2104427 (86.07%)
# reads that failed to align: 191414 (7.83%)
# reads with alignments suppressed due to -m: 149057 (6.10%)
# reads processed: 13231701                                                                                                                                     
# 1320                                                               
# reads with at least one reported alignment: 11548572 (87.28%)                                                                                                                                                                           
# reads that failed to align: 1131652 (8.55%)   
# reads with alignments suppressed due to -m: 551477 (4.17%)
# Reported 11548572 alignments to 1 output stream(s)    
