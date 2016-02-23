#!/bin/sh



ecoli=/p/keles/ChIPexo/volume6/K12/pther_methods/gem/ecoli.size
gemdir=$HOME/Desktop/Docs/Code/lib/gem
gem=$gemdir/gem.jar
readdistr=$gemdir/Read_Distribution_ChIP-exo.txt

basedir=/p/keles/ChIPexo/volume7/Landick
genome=$basedir/index/E.coli_K_12_MG1655_GENOME.fa

k=6
nCore=20
size=4400000

FDR=FDR10

beddir=/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo/bedfiles
outdir=/p/keles/ChIPexo/volume6/K12/other_methods/gem
peakdir=$outdir/peaks

# echo $exofile
# echo $outfile
# echo $peakfile

 # --g $ecoli

exofile=$beddir/edsn1311_Sig70.sort.bed
outfile=$outdir/edsn1311_Sig70_$FDR
peakfile=$peakdir/edsn1311_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED

exofile=$beddir/edsn1314_Sig70.sort.bed
outfile=$outdir/edsn1314_Sig70_$FDR
peakfile=$peakdir/edsn1314_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED

exofile=$beddir/edsn1317_Sig70.sort.bed
outfile=$outdir/edsn1317_Sig70_$FDR
peakfile=$peakdir/edsn1317_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED

exofile=$beddir/edsn1320_Sig70.sort.bed
outfile=$outdir/edsn1320_Sig70_$FDR
peakfile=$peakdir/edsn1320_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED

exofile=$beddir/edsn931_Sig70.sort.bed
outfile=$outdir/edsn931_Sig70_$FDR
peakfile=$peakdir/edsn931_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED

exofile=$beddir/edsn933_Sig70.sort.bed
outfile=$outdir/edsn933_Sig70_$FDR
peakfile=$peakdir/edsn933_Sig70_peaks_$FDR.txt

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED
