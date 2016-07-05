#!/bin/sh

sizefile=/p/keles/ChIPexo/volume7/Landick/K12/K12_size
outdir=data/SCC_curves
outdir2=data/ChIPseq_QC
maxShift=300



# landick aero sig70 - SE ChIP-seq
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET/aerobic_vs_anaerobic_bed
isPET=FALSE

rscripts/scripts/SCC_curve_and_QC_ChIPseq_SE.R $dr/run80.sigma70+O2_I_PA.s_6_sequence.eland.Umatch.txt.bed $outdir/landick_sig70+O2_SE_SCC.txt $outdir2/landick_sig70+O2_SE_QC.txt $sizefile $isPET $maxShift

rscripts/scripts/SCC_curve_and_QC_ChIPseq_SE.R $dr/run80.sigma70-O2_IP_A.s_6_sequence.txt.eland.Umatch.txt.bed $outdir/landick_sig70-O2_SE_SCC.txt $outdir2/landick_sig70-O2_SE_QC.txt $sizefile $isPET $maxShift


# landick aero sig70 - PE ChIP-seq
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_PET/aerobic_vs_anaerobic_bed
isPET=TRUE

rscripts/scripts/SCC_curve_and_QC_ChIPseq_PE.R $dr/run208.sigma70-minuso2a_TGACCA_L003.paired_Umatch_eland.tab.bed $dr/run208.sigma70-minuso2a_TGACCA_L003.paired_Umatch_eland.tab_coord.txt $outdir/landick_sig70-O2_PE_SCC.txt $outdir2/landick_sig70-O2_PE_QC.txt $sizefile $isPET $maxShift

rscripts/scripts/SCC_curve_and_QC_ChIPseq_PE.R $dr/run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab.bed $dr/run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab_coord.txt $outdir/landick_sig70+O2_PE_SCC.txt $outdir2/landick_sig70+O2_PE_QC.txt $sizefile $isPET $maxShift
