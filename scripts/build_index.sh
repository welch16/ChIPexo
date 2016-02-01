#!/bin/sh

# script that builds the index to use for Bob Landick's ChIP-exo data

dr=/p/keles/ChIPexo/volume7/Landick/index
file=RL3000_rl111612jg.fas

bowtie-build -f $dr/$file $dr/RL3000
