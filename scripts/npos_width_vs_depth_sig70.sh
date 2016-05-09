#!/bin/sh

N=10000
M=1000
indir=data/ChIPexo_QC_runs_sig70
outdir=data/npos_width_reg_sig70

for f in "$indir"/*RData
do
  g=${f/$indir/$outdir}
  g=${g/.RData/_factors.RData}
  rscripts/scripts/npos_width_vs_depth_regression.R $f $g $M $N
done

