#!/bin/sh

outdir=data/SCC_curves
maxShift=300
isPET=FALSE

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

rscripts/scripts/SCC_curve.R $dr/ERR336942.sort.bam $outdir/carroll_mouse_Rep1_SCC.txt mm9 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ERR336956.sort.bam $outdir/carroll_mouse_Rep2_SCC.txt mm9 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ERR336935.sort.bam $outdir/carroll_mouse_Rep3_SCC.txt mm9 $isPET $maxShift

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

rscripts/scripts/SCC_curve.R $dr/ERR336933.sort.bam $outdir/carroll_human_Rep1_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ERR336950.sort.bam $outdir/carroll_human_Rep2_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ERR336938.sort.bam $outdir/carroll_human_Rep3_SCC.txt hg19 $isPET $maxShift

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data

rscripts/scripts/SCC_curve.R $dr/CTCF.bam $outdir/pugh_CTCF_SCC.txt hg19 $isPET $maxShift

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/SCC_curve.R $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam $outdir/chipnexus_embryo_dorsal_Rep1_SCC.txt dm3 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam $outdir/chipnexus_embryo_dorsal_Rep2_SCC.txt dm3 $isPET $maxShift

rscripts/scripts/SCC_curve.R $dr/ChIPnexus_embryo_twist_rep1.sort.bam $outdir/chipnexus_embryo_twist_Rep1_SCC.txt dm3 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ChIPnexus_embryo_twist_rep2.sort.bam $outdir/chipnexus_embryo_twist_Rep2_SCC.txt dm3 $isPET $maxShift

rscripts/scripts/SCC_curve.R $dr/ChIPnexus_K562_TBP_rep1.sort.bam $outdir/chipnexus_K562_TBP_Rep1_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ChIPnexus_K562_TBP_rep2.sort.bam $outdir/chipnexus_K562_TBP_Rep2_SCC.txt hg19 $isPET $maxShift

rscripts/scripts/SCC_curve.R $dr/ChIPnexus_S2_Max_rep1.sort.bam $outdir/chipnexus_S2_Max_Rep1_SCC.txt dm3 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ChIPnexus_S2_Max_rep2.sort.bam $outdir/chipnexus_S2_Max_Rep2_SCC.txt dm3 $isPET $maxShift

rscripts/scripts/SCC_curve.R $dr/ChIPnexus_S2_MyC_rep1.sort.bam $outdir/chipnexus_S2_MyC_Rep1_SCC.txt dm3 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/ChIPnexus_S2_MyC_rep2.sort.bam $outdir/chipnexus_S2_MyC_Rep2_SCC.txt dm3 $isPET $maxShift

# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/SCC_curve.R $dr/TBP_K562_Rep1.sort.bam $outdir/venters_TBP_K562_Rep1_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/TBP_K562_Rep2.sort.bam $outdir/venters_TBP_K562_Rep2_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/TBP_K562_Rep3.sort.bam $outdir/venters_TBP_K562_Rep3_SCC.txt hg19 $isPET $maxShift

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

rscripts/scripts/SCC_curve.R $dr/IMR90_GR_chip-exo.sort.bam $outdir/meijsing_GR_IMR90_Rep1_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/K562_GR_chip-exo.sort.bam $outdir/meijsing_GR_K562_Rep1_SCC.txt hg19 $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/U2OS_GR_chip-exo.sort.bam $outdir/meijsing_GR_U2OS_Rep1_SCC.txt hg19 $isPET $maxShift


# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic
sizefile=/p/keles/ChIPexo/volume7/Landick/K12/K12_size

rscripts/scripts/SCC_curve.R $dr/edsn931_Sig70.sort.bam $outdir/landick_sig70_931_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn933_Sig70.sort.bam $outdir/landick_sig70_933_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn935_Sig70.sort.bam $outdir/landick_sig70_935_SCC.txt $sizefile $isPET $maxShift 
rscripts/scripts/SCC_curve.R $dr/edsn937_Sig70.sort.bam $outdir/landick_sig70_937_SCC.txt $sizefile $isPET $maxShift

# landick rif sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

rscripts/scripts/SCC_curve.R $dr/edsn1311_Sig70.sort.bam $outdir/landick_sig70_1311_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1314_Sig70.sort.bam $outdir/landick_sig70_1314_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1317_Sig70.sort.bam $outdir/landick_sig70_1317_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1320_Sig70.sort.bam $outdir/landick_sig70_1320_SCC.txt $sizefile $isPET $maxShift

# landick rif sig70 - SE ChIP-seq
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET/rif_treatment

rscripts/scripts/SCC_curve.R $dr/edsn1396_Sig70.sort.bam $outdir/landick_sig70_1396_SE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1398_Sig70.sort.bam $outdir/landick_sig70_1398_SE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1400_Sig70.sort.bam $outdir/landick_sig70_1400_SE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1402_Sig70.sort.bam $outdir/landick_sig70_1402_SE_SCC.txt $sizefile $isPET $maxShift

# landick rif sig70 - PE ChIP-seq
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_PET/rif_treatment
isPET=TRUE

rscripts/scripts/SCC_curve.R $dr/edsn1396_Sig70.sort.bam $outdir/landick_sig70_1396_PE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1398_Sig70.sort.bam $outdir/landick_sig70_1398_PE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1400_Sig70.sort.bam $outdir/landick_sig70_1400_PE_SCC.txt $sizefile $isPET $maxShift
rscripts/scripts/SCC_curve.R $dr/edsn1402_Sig70.sort.bam $outdir/landick_sig70_1402_PE_SCC.txt $sizefile $isPET $maxShift




