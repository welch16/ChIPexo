#!/bin/sh


ecoli=inst/gem_analysis/ecoli.size
gemdir=$HOME/Desktop/Docs/Code/lib/gem
basedir=/p/keles/ChIPexo/volume7/Landick
gem=$gemdir/gem.jar
readdistr=$gemdir/Read_Distribution_ChIP-exo.txt
genome=$basedir/index/RL3000_rl111612.fa
k=6
nCore=20
size=4400000

exofile=$basedir/ChIPexo/bedfiles/edsn1311_Sig70.sort.bed
outfile=inst/gem_analysis/edsn1311_Sig70

 # --genome $genome

java -Xmx10G -jar $gem --t $nCore --g $ecoli --d $readdistr --s $size --exptX $exofile --k $k --out $outfile

exofile=$basedir/ChIPexo/bedfiles/edsn1314_Sig70.sort.bed
outfile=inst/gem_analysis/edsn1314_Sig70

java -Xmx10G -jar $gem --t $nCore --g $ecoli --d $readdistr --s $size --exptX $exofile --k $k --out $outfile

exofile=$basedir/ChIPexo/bedfiles/edsn1317_Sig70.sort.bed
outfile=inst/gem_analysis/edsn1317_Sig70

java -Xmx10G -jar $gem --t $nCore --g $ecoli --d $readdistr --s $size --exptX $exofile --k $k --out $outfile

exofile=$basedir/ChIPexo/bedfiles/edsn1320_Sig70.sort.bed
outfile=inst/gem_analysis/edsn1320_Sig70

java -Xmx10G -jar $gem --t $nCore --g $ecoli --d $readdistr --s $size --exptX $exofile --k $k --out $outfile
