#!/bin/sh


index=/p/keles/ChIPexo/volume7/Landick/index/E.Coli_K-12
mc=24
maxAlign=1
minins=100
maxins=500

basedir=/p/keles/ChIPexo
indir=$basedir/volume8/RawData_from_Biotech_04Dec15
indir=$indir/run423
outdir=$basedir/volume7
outdir=$outdir/Landick/K12/ChIPseq_PET/rif_treatment

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



# Aligning ChIP seq PET files:
# Fragment length alowed in the 100 - 500 range
# Index file: /p/keles/ChIPexo/volume7/Landick/index/E.Coli_K-12
# Threads used: 24
# Allowing reads to align to at most 1 positions
# # reads processed: 11825205
# # reads with at least one reported alignment: 6722511 (56.85%)
# # reads that failed to align: 4877126 (41.24%)
# # reads with alignments suppressed due to -m: 225568 (1.91%)
# Reported 6722511 paired-end alignments to 1 output stream(s)
# # reads processed: 14985772
# # reads with at least one reported alignment: 8269460 (55.18%)
# # reads that failed to align: 6482333 (43.26%)
# # reads with alignments suppressed due to -m: 233979 (1.56%)
# Reported 8269460 paired-end alignments to 1 output stream(s)
# # reads processed: 9005641
# # reads with at least one reported alignment: 5821361 (64.64%)
# # reads that failed to align: 3011056 (33.44%)
# # reads with alignments suppressed due to -m: 173224 (1.92%)
# Reported 5821361 paired-end alignments to 1 output stream(s)
# # reads processed: 14002998
# # reads with at least one reported alignment: 8427013 (60.18%)
# # reads that failed to align: 5348751 (38.20%)
# # reads with alignments suppressed due to -m: 227234 (1.62%)
# Reported 8427013 paired-end alignments to 1 output stream(s)
# # reads processed: 16923160
# # reads with at least one reported alignment: 10374674 (61.30%)
# # reads that failed to align: 6320137 (37.35%)
# # reads with alignments suppressed due to -m: 228349 (1.35%)
# Reported 10374674 paired-end alignments to 1 output stream(s)
