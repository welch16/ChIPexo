#!/bin/sh

indir=data/ChIPexo_QC_runs
outdir=figs/ARC_factors

rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_human_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_human_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_human_Rep3.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_mouse_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_mouse_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/carroll_mouse_Rep3.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_embryo_dorsal_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_embryo_dorsal_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_embryo_twist_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_embryo_twist_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_K562_TBP_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_K562_TBP_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_S2_Max_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_S2_Max_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_S2_MyC_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/chipnexus_S2_MyC_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_1311.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_1314.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_1317.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_1320.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_931.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_933.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_935.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/landick_sig70_937.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/meijsing_GR_IMR90_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/meijsing_GR_K562_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/meijsing_GR_U2OS_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/pugh_CTCF.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/venters_TBP_K562_Rep1.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/venters_TBP_K562_Rep2.RData $outdir
rscripts/scripts/arc_vs_urcr_factor_test.R $indir/venters_TBP_K562_Rep3.RData $outdir
