#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/meme
cores=20
minw=20
maxw=35
ms=1000000
nmot=5

base=$PWD

cd $indir

meme edsn1396_Sig70_sequences.s -oc edsn1396_meme_results_anr -dna -mod anr -p $cores -minw $minw -maxw $maxw -revcomp -maxsize $ms -nmotifs $nmot

meme edsn1398_Sig70_sequences.s -oc edsn1398_meme_results_anr -dna -mod anr -p $cores -minw $minw -maxw $maxw -revcomp -maxsize $ms -nmotifs $nmot

meme edsn1400_Sig70_sequences.s -oc edsn1400_meme_results_anr -dna -mod anr -p $cores -minw $minw -maxw $maxw -revcomp -maxsize $ms -nmotifs $nmot

meme edsn1402_Sig70_sequences.s -oc edsn1402_meme_results_anr -dna -mod anr -p $cores -minw $minw -maxw $maxw -revcomp -maxsize $ms -nmotifs $nmot


cd $base
