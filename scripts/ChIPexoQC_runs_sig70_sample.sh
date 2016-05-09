#!/bin/sh

outdir=data/ChIPexo_QC_runs_sig70
seed=$1

basedir=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo
dir1=$basedir/aerobic_vs_anaerobic
dir2=$basedir/rif_treatment


for i in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
do
    for f in $dir1/*sort.bam
    do
	g=${f/$dir1/$outdir}
	g=${g/.sort.bam/_subsample$i.RData}
	rscripts/scripts/ChIPexoQual_run.R $f $g $i $seed
    done
done

for i in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
do
    for f in $dir2/*sort.bam
    do
	g=${f/$dir2/$outdir}
	g=${g/.sort.bam/_subsample$i.RData}
	rscripts/scripts/ChIPexoQual_run.R $f $g $i $seed
    done
done


