#!/bin/sh

pipeline=./rscripts/figures/fig5/fig5_ChIPexoQual_pipeline.R
outdr=./data/figures/fig5/ecoli/1000
cores=4
nregions=800
ntimes=1000
basedr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo

#########################################################

# Landick E1
indr=$basedr/aerobic_vs_anaerobic

$pipeline --readsfile $indr/edsn931_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_931_stats.tsv \
	  --scoresfile $outdr/Landick_931_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn933_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_933_stats.tsv \
	  --scoresfile $outdr/Landick_933_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn935_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_935_stats.tsv \
	  --scoresfile $outdr/Landick_935_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn937_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_937_stats.tsv \
	  --scoresfile $outdr/Landick_937_scores.tsv \
	  --cores $cores

# Landick E2
indr=$basedr/rif_treatment

$pipeline --readsfile $indr/edsn1311_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_1311_stats.tsv \
	  --scoresfile $outdr/Landick_1311_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn1314_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_1314_stats.tsv \
	  --scoresfile $outdr/Landick_1314_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn1317_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_1317_stats.tsv \
	  --scoresfile $outdr/Landick_1317_scores.tsv \
	  --cores $cores
$pipeline --readsfile $indr/edsn1320_Sig70.sort.bam \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/Landick_1320_stats.tsv \
	  --scoresfile $outdr/Landick_1320_scores.tsv \
	  --cores $cores
