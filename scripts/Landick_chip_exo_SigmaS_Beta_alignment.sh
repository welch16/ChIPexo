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

bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq $outdir1/edsn931_Sig70.sam

# bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1003b-No-Ab-Lib_TGACCA_L007_R1.fastq

bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1010-Sig-S-Lib_ACTTGA_L007_R1.fastq $outdir1/edsn938_SigS.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1003b-Sig-S-Lib_CTTGTA_L007_R1.fastq $outdir1/edsn936_SigS.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1003a-Sig-S-Lib_TAGCTT_L007_R1.fastq $outdir1/edsn934_SigS.sam
bowtie -q -m 1 -l 55 -k 1 -5 3 -3 40 --best -S -p $mc $index $indir1/run299.456-166-1002-Sig-S-Lib_GCCAAT_L007_R1.fastq $outdir1/edsn932_SigS.sam

## new files

echo --- new files ---

# bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1311-cult-1197-ChIPexo-0minRif-Sigma70_CGATGT_L001_R1.fastq $outdir2/edsn1311_Sig70.sam


bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1321-cult-1202-ChIPexo-20minRif-BetaPrimeFLAG_CTTGTA_L001_R1.fastq $outdir2/edsn1321_BetaPF.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1319-cult-1202-ChIPexo-20minRif-Beta_TAGCTT_L001_R1.fastq $outdir2/edsn1319_Beta.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1318-cult-1202-ChIPexo-0minRif-BetaPrimeFLAG_GATCAG_L001_R1.fastq $outdir2/edsn1318_BetaPF.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1316-cult-1202-ChIPexo-0minRif-Beta_CAGATC_L001_R1.fastq $outdir2/edsn1316_Beta.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1315-cult-1197-ChIPexo-20minRif-BetaPrimeFLAG_GCCAAT_L001_R1.fastq $outdir2/edsn1315_BetaPF.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1313-cult-1197-ChIPexo-20minRif-Beta_TGACCA_L001_R1.fastq $outdir2/edsn1313_Beta.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1312-cult-1197-ChIPexo-0minRif-BetaPrimeFLAG_TTAGGC_L001_R1.fastq $outdir2/edsn1312_BetaPF.sam
bowtie -q -m 1 -v 2 --best -p $mc -S $index $indir2/run423.edsn-1310-cult-1197-ChIPexo-0minRif-Beta_ATCACG_L001_R1.fastq $outdir2/edsn1310_Beta.sam
