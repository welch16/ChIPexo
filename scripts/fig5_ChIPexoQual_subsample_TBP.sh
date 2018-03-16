#!/bin/sh

pipeline=./rscripts/figures/fig5/fig5_ChIPexoQual_subsampling.R
outdr=./data/figures/fig5
cores=22
nregions=1000
ntimes=1000
basedr=/p/keles/ChIPexo/volume4

#########################################################

indr=$basedr/venters_data/sortbam

$pipeline --readsfile $indr/TBP_K562_Rep1.sort.bam \
	  --depths 20,30,40,50 --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/venters_TBP_K562_Rep1_stats.tsv \
	  --scoresfile $outdr/venters_TBP_K562_Rep1_scores.tsv \
	  --cores $cores

$pipeline --readsfile $indr/TBP_K562_Rep2.sort.bam \
	  --depths 20,30,40,50 --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/venters_TBP_K562_Rep2_stats.tsv \
	  --scoresfile $outdr/venters_TBP_K562_Rep2_scores.tsv \
	  --cores $cores

$pipeline --readsfile $indr/TBP_K562_Rep3.sort.bam \
	  --depths 20,30,40,50 --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/venters_TBP_K562_Rep3_stats.tsv \
	  --scoresfile $outdr/venters_TBP_K562_Rep3_scores.tsv \
	  --cores $cores

indr=$basedr/zeitlinger_data/bam/sortbam


$pipeline --readsfile $indr/ChIPnexus_K562_TBP_rep2.sort.bam \
	  --depths 20,30,40,50 --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/zeitlinger_TBP_K562_Rep2_stats.tsv \
	  --scoresfile $outdr/zeitlinger_TBP_K562_Rep2_scores.tsv \
	  --cores $cores
