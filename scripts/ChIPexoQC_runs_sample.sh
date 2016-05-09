#!/bin/sh

outdir=data/ChIPexo_QC_runs
N1=20000000
N2=50000000
seed=$1

## venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1_subsample20.RData $N1 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2_subsample20.RData $N1 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3_subsample20.RData $N1 $seed


rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1_subsample50.RData $N2 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2_subsample50.RData $N2 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3_subsample50.RData $N2 $seed

N3=30000000
N4=40000000

rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1_subsample30.RData $N3 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2_subsample30.RData $N3 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3_subsample30.RData $N3 $seed

rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1_subsample40.RData $N4 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2_subsample40.RData $N4 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3_subsample40.RData $N4 $seed


## chip nexus K562

dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2_subsample20.RData $N1 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2_subsample50.RData $N2 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2_subsample30.RData $N3 $seed
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2_subsample40.RData $N4 $seed
