#!/bin/sh

index=$HOME/Desktop/Docs/Code/lib/bowtie2-2.2.6/indexes/ecoli/NC_000913
mc=24
maxAlign=1
minins=100
maxins=500

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir=$indir/run423
outdir=$basedir/volume7
outdir=$outdir/rif_treatment

echo Aligning ChIP seq PET files:
echo Fragment length alowed in the $minins - $maxins range
echo Index file: $index
echo Threads used: $mc
echo Allowing reads to align to at most $maxAlign positions

# run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R1.fastq.gz
# run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R2.fastq.gz

bowtie2 --threads $mc -k $maxAlign -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R1.fastq.gz -2 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R2.fastq.gz -S $outdir/edsn1396_Sig70.sam

# run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R1.fastq.gz
# run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R2.fastq.gz

bowtie2 --threads $mc -k $maxAlign -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R1.fastq.gz -2 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R2.fastq.gz -S $outdir/edsn1398_Sig70.sam

# run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R1.fastq.gz
# run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R2.fastq.gz

bowtie2 --threads $mc -k $maxAlign -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R1.fastq.zg -2 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R2.fastq.gz -S $outdir/edsn1400_Sig70.sam

# run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R1.fastq.gz
# run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R2.fastq.gz

bowtie2 --threads $mc -k $maxAlign -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R1.fastq.gz -2 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R2.fastq.gz -S $outdir/edsn1402_Sig70.sam

# run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R1.fastq.gz
# run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R2.fastq.gz

bowtie2 --threads $mc -k $maxAlign -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R1.fastq.gz -2 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R2.fastq.gz -S $outdir/edsn1369_Input.sam
