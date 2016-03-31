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

beddir=/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo/bedfiles
outdir=/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test

# echo $exofile
# echo $outfile
# echo $peakfile

 # --g $ecoli

exofile=$beddir/edsn931_Sig70.sort.bed
outfile=$outdir/edsn931_Sig70_$FDR
peakfile=/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test/peak_to_check_gem.txt
outfile=$outdir/gem_example

java -Xmx10G -jar $gem --t $nCore --d $readdistr --s $size --subf $peakfile --exptX $exofile --k $k --out $outfile --outBED
