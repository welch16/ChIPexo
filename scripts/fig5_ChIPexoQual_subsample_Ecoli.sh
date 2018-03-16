pipeline=./rscripts/figures/fig5/fig5_ChIPexoQual_subsampling.R
cores=22
nregions=800
ntimes=1000
outdr=./data/figures/fig5/ecoli/"$nregions"
basedr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo

#########################################################


# Landick E1
indr=$basedr/aerobic_vs_anaerobic

$pipeline --readsfile $indr/edsn931_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn931_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn931_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn933_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn933_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn933_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn935_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn935_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn935_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn937_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn937_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn937_Sig70_subsample_scores.tsv --cores $cores


# Landick E2
indr=$basedr/rif_treatment

$pipeline --readsfile $indr/edsn1311_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn1311_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn1311_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn1314_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn1314_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn1314_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn1317_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn1317_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn1317_Sig70_subsample_scores.tsv --cores $cores
$pipeline --readsfile $indr/edsn1320_Sig70.sort.bam \
	  --depths 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1 \
	  --nregions $nregions --ntimes $ntimes \
	  --statsfile $outdr/edsn1320_Sig70_subsample_stats.tsv \
	  --scoresfile $outdr/edsn1320_Sig70_subsample_scores.tsv --cores $cores


