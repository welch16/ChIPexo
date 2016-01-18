#!/bin/sh


index=/p/keles/ChIPexo/volume7/Landick/index/RL3000
mc=24
maxAlign=1
minins=100
maxins=500

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir=$indir/run423
outdir=$basedir/volume7
outdir=$outdir/Landick/ChIPseq_PET/rif_treatment

echo Aligning ChIP seq PET files:
echo Fragment length alowed in the $minins - $maxins range
echo Index file: $index
echo Threads used: $mc
echo Allowing reads to align to at most $maxAlign positions

# run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R1.fastq.gz
# run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R2.fastq.gz

bowtie -q -v 2 -a -m 1 -p $mc -S -I $minins -X $maxins $index -1 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R1.fastq -2 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R2.fastq $outdir/edsn1396_Sig70.sam

# bowtie2 --threads $mc -k $maxAlign --ff -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R1.fastq.gz -2 $indir/run423.edsn-1396-cult-1197-ChIPseqIP-0minRif-Sigma70_GTGAAA_L002_R2.fastq.gz -S $outdir/edsn1396_Sig70.sam

# run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R1.fastq.gz
# run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R2.fastq.gz

bowtie -q -v 2 -a -m 1 -p $mc -S -I $minins -X $maxins $index -1 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R1.fastq -2 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R2.fastq -S $outdir/edsn1398_Sig70.sam

# bowtie2 --threads $mc -k $maxAlign --ff -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R1.fastq.gz -2 $indir/run423.edsn-1398-cult-1197-ChIPseqIP-20minRif-Sigma70_GTTTCG_L002_R2.fastq.gz -S $outdir/edsn1398_Sig70.sam

# run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R1.fastq.gz
# run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R2.fastq.gz

bowtie -q -v 2 -a -m 1 -p $mc -S -I $minins -X $maxins $index -1 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R1.fastq -2 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R2.fastq -S $outdir/edsn1400_Sig70.sam


# bowtie2 --threads $mc -k $maxAlign --ff -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R1.fastq.gz -2 $indir/run423.edsn-1400-cult-1202-ChIPseqIP-0minRif-Sigma70_GAGTGG_L002_R2.fastq.gz -S $outdir/edsn1400_Sig70.sam

# run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R1.fastq.gz
# run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R2.fastq.gz

bowtie -q -v 2 -a -m 1 -p $mc -S -I $minins -X $maxins $index -1 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R1.fastq -2 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R2.fastq -S $outdir/edsn1402_Sig70.sam

# bowtie2 --threads $mc -k $maxAlign --ff -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R1.fastq.gz -2 $indir/run423.edsn-1402-cult-1202-ChIPseqIP-20minRif-Sigma70_ATTCCT_L002_R2.fastq.gz -S $outdir/edsn1402_Sig70.sam

# run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R1.fastq.gz
# run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R2.fastq.gz

bowtie -q -v 2 -a -m 1 -p $mc -S -I $minins -X $maxins $index -1 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R1.fastq -2 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R2.fastq -S $outdir/edsn1369_Input.sam

# bowtie2 --threads $mc -k $maxAlign --ff -I $minins -X $maxins -x $index -1 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R1.fastq.gz -2 $indir/run423.edsn-1369-cult-1204-ChIPseqInput-midlog-Beta_TTAGGC_L002_R2.fastq.gz -S $outdir/edsn1369_Input.sam


# bowtie -q -v 2 -a -m 1 -p $mc -S $index $indir2/run423.edsn-1320-cult-1202-ChIPexo-20minRif-Sigma70_GGCTAC_L001_R1.fastq $outdir2/edsn1320_Sig70.sam
