#!/bin/sh

basedir=/p/keles/ChIPexo/volume6/K12
base=$PWD

meme=./meme_sig70/sig70_motif_oops_test/meme.txt


cd $basedir

for e in 1311 1314 1317 1320 931 933 935 937
do
    seq=./meme/edsn$e"_Sig70_sequences.s"
    out=./meme_sig70/edsn$e"_fimo_motif_all_test"
    fimo -oc $out $meme $seq
done

cd $base


    # for p in 1396 1398 1400 1402
    # do
    # 	out=exo_fimo/edsn$e"_motif_"$p
    # 	meme="edsn"$p"_meme_results/meme.txt"
    # 	seq="edsn"$e"_Sig70_sequences.s"
    # 	# echo "---"
    # 	# echo $out
    # 	# echo $meme
    # 	# echo $seq
    # 	fimo -oc $out --motif 2 $meme $seq
    # done
