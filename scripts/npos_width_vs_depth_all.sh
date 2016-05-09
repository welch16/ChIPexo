#!/bin/sh

N=10000
M=1000
indir=data/ChIPexo_QC_runs
outdir=data/npos_width_reg

for f in "$indir"/*RData
do
  g=${f/$indir/$outdir}
  g=${g/.RData/_factors.RData}
  rscripts/scripts/npos_width_vs_depth_regression.R $f $g $M $N
done

