#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/meme

base=$PWD

cd $indir

for e in 1311 1314 1317 1320 931 933 935 937
do
    for p in 1396 1398 1400 1402
    do
	out=exo_fimo/edsn$e"_motif_"$p
	meme="edsn"$p"_meme_results/meme.txt"
	seq="edsn"$e"_Sig70_sequences.s"
	# echo "---"
	# echo $out
	# echo $meme
	# echo $seq
	fimo -oc $out --motif 2 $meme $seq
    done
done

cd $base
