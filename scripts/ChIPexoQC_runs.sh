#!/bin/sh

outdir=data/ChIPexo_QC_runs

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

rscripts/scripts/ChIPexoQual_run.R $dr/ERR336942.sort.bam $outdir/carroll_mouse_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ERR336956.sort.bam $outdir/carroll_mouse_Rep2.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ERR336935.sort.bam $outdir/carroll_mouse_Rep3.RData 

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

rscripts/scripts/ChIPexoQual_run.R $dr/ERR336933.sort.bam $outdir/carroll_human_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ERR336950.sort.bam $outdir/carroll_human_Rep2.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ERR336938.sort.bam $outdir/carroll_human_Rep3.RData 

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data

rscripts/scripts/ChIPexoQual_run.R $dr/CTCF.bam $outdir/pugh_CTCF.RData 

# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic

rscripts/scripts/ChIPexoQual_run.R $dr/edsn931_Sig70.sort.bam $outdir/landick_sig70_931.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn933_Sig70.sort.bam $outdir/landick_sig70_933.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn935_Sig70.sort.bam $outdir/landick_sig70_935.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn937_Sig70.sort.bam $outdir/landick_sig70_937.RData 

# landick rif sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

rscripts/scripts/ChIPexoQual_run.R $dr/edsn1311_Sig70.sort.bam $outdir/landick_sig70_1311.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn1314_Sig70.sort.bam $outdir/landick_sig70_1314.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn1317_Sig70.sort.bam $outdir/landick_sig70_1317.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/edsn1320_Sig70.sort.bam $outdir/landick_sig70_1320.RData 

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam $outdir/chipnexus_embryo_dorsal_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam $outdir/chipnexus_embryo_dorsal_Rep2.RData 

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_embryo_twist_rep1.sort.bam $outdir/chipnexus_embryo_twist_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_embryo_twist_rep2.sort.bam $outdir/chipnexus_embryo_twist_Rep2.RData 

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep1.sort.bam $outdir/chipnexus_K562_TBP_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2.RData 

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_S2_Max_rep1.sort.bam $outdir/chipnexus_S2_Max_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_S2_Max_rep2.sort.bam $outdir/chipnexus_S2_Max_Rep2.RData 

rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_S2_MyC_rep1.sort.bam $outdir/chipnexus_S2_MyC_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/ChIPnexus_S2_MyC_rep2.sort.bam $outdir/chipnexus_S2_MyC_Rep2.RData 

# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3.RData 

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

rscripts/scripts/ChIPexoQual_run.R $dr/IMR90_GR_chip-exo.sort.bam $outdir/meijsing_GR_IMR90_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/K562_GR_chip-exo.sort.bam $outdir/meijsing_GR_K562_Rep1.RData 
rscripts/scripts/ChIPexoQual_run.R $dr/U2OS_GR_chip-exo.sort.bam $outdir/meijsing_GR_U2OS_Rep1.RData 



