#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/meme
cores=20
minw=30
maxw=40

base=$PWD

cd $indir

mkdir edsn1396_meme_results
meme edsn1396_Sig70_sequences.s -oc edsn1396_meme_results -text -dna -mod zoops -p $cores -minw $minw -maxw $maxw -revcomp -maxsize 200000

mkdir edsn1398_meme_results
meme edsn1398_Sig70_sequences.s -oc edsn1398_meme_results -text -dna -mod zoops -p $cores -minw $minw -maxw $maxw -revcomp -maxsize 200000

mkdir edsn1400_meme_results
meme edsn1400_Sig70_sequences.s -oc edsn1400_meme_results -text -dna -mod zoops -p $cores -minw $minw -maxw $maxw -revcomp -maxsize 200000

mkdir edsn1402_meme_results
meme edsn1402_Sig70_sequences.s -oc edsn1402_meme_results -text -dna -mod zoops -p $cores -minw $minw -maxw $maxw -revcomp -maxsize 200000

cd $base
