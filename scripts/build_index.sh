#!/bin/sh

# script that builds the index to use for Bob Landick's ChIP-exo data

dr=/p/keles/ChIPexo/volume7/Landick/index
file=E.coli_K-12_MG1655_GENOME.fa

bowtie-build -f $dr/$file $dr/E.Coli_K-12
