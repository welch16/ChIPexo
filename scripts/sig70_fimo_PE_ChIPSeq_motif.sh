#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/meme

base=$PWD

cd $indir

fimo -oc edsn1396_fimo_results --motif 2 edsn1396_meme_results/meme.txt edsn1396_Sig70_sequences.s

fimo -oc edsn1398_fimo_results --motif 2 edsn1398_meme_results/meme.txt edsn1398_Sig70_sequences.s

fimo -oc edsn1400_fimo_results --motif 2 edsn1400_meme_results/meme.txt edsn1400_Sig70_sequences.s

fimo -oc edsn1402_fimo_results --motif 2 edsn1402_meme_results/meme.txt edsn1402_Sig70_sequences.s


cd $base
