#!/bin/sh

outdir=data/SCC_spp_curves
maxShift=300
shift=1,$maxShift
cores=12

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336942.sort.bam --outfile $outdir/carroll_mouse_Rep1_SCC.txt --scc.range $shift --ncluster 12 
rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336956.sort.bam --outfile $outdir/carroll_mouse_Rep2_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336935.sort.bam --outfile $outdir/carroll_mouse_Rep3_SCC.txt --scc.range $shift --ncluster 12

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336933.sort.bam --outfile $outdir/carroll_human_Rep1_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336950.sort.bam --outfile $outdir/carroll_human_Rep2_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ERR336938.sort.bam --outfile $outdir/carroll_human_Rep3_SCC.txt --scc.range $shift --ncluster 12

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data

rscripts/scripts/spp_QC_run.R --chipfile $dr/CTCF.bam --outfile $outdir/pugh_CTCF_SCC.txt --scc.range $shift --ncluster 12

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam --outfile $outdir/chipnexus_embryo_dorsal_Rep1_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam --outfile $outdir/chipnexus_embryo_dorsal_Rep2_SCC.txt --scc.range $shift --ncluster 12 

rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_embryo_twist_rep1.sort.bam --outfile $outdir/chipnexus_embryo_twist_Rep1_SCC.txt  --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_embryo_twist_rep2.sort.bam --outfile $outdir/chipnexus_embryo_twist_Rep2_SCC.txt  --scc.range $shift --ncluster 12


rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_K562_TBP_rep1.sort.bam --outfile $outdir/chipnexus_K562_TBP_Rep1_SCC.txt  --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_K562_TBP_rep2.sort.bam --outfile $outdir/chipnexus_K562_TBP_Rep2_SCC.txt  --scc.range $shift --ncluster 12


rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_S2_Max_rep1.sort.bam --outfile $outdir/chipnexus_S2_Max_Rep1_SCC.txt  --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_S2_Max_rep2.sort.bam --outfile $outdir/chipnexus_S2_Max_Rep2_SCC.txt  --scc.range $shift --ncluster 12


rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_S2_MyC_rep1.sort.bam --outfile $outdir/chipnexus_S2_MyC_Rep1_SCC.txt  --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/ChIPnexus_S2_MyC_rep2.sort.bam --outfile $outdir/chipnexus_S2_MyC_Rep2_SCC.txt  --scc.range $shift --ncluster 12


# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/spp_QC_run.R --chipfile $dr/TBP_K562_Rep1.sort.bam --outfile $outdir/venters_TBP_K562_Rep1_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/TBP_K562_Rep2.sort.bam --outfile $outdir/venters_TBP_K562_Rep2_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/TBP_K562_Rep3.sort.bam --outfile $outdir/venters_TBP_K562_Rep3_SCC.txt --scc.range $shift --ncluster 12

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

rscripts/scripts/spp_QC_run.R --chipfile $dr/IMR90_GR_chip-exo.sort.bam --outfile $outdir/meijsing_GR_IMR90_Rep1_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/K562_GR_chip-exo.sort.bam --outfile $outdir/meijsing_GR_K562_Rep1_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/U2OS_GR_chip-exo.sort.bam --outfile $outdir/meijsing_GR_U2OS_Rep1_SCC.txt --scc.range $shift --ncluster 12


# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic

rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn931_Sig70.sort.bam --outfile $outdir/landick_sig70_931_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn933_Sig70.sort.bam --outfile $outdir/landick_sig70_933_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn935_Sig70.sort.bam --outfile $outdir/landick_sig70_935_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn937_Sig70.sort.bam --outfile $outdir/landick_sig70_937_SCC.txt --scc.range $shift --ncluster 12

# landick rif sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1311_Sig70.sort.bam --outfile $outdir/landick_sig70_1311_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1314_Sig70.sort.bam --outfile $outdir/landick_sig70_1314_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1317_Sig70.sort.bam --outfile $outdir/landick_sig70_1317_SCC.txt --scc.range $shift --ncluster 12
rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1320_Sig70.sort.bam --outfile $outdir/landick_sig70_1320_SCC.txt --scc.range $shift --ncluster 12

# # landick rif sig70 - SE ChIP-seq
# dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET/rif_treatment

# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1396_Sig70.sort.bam --outfile $outdir/landick_sig70_1396_SE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1398_Sig70.sort.bam --outfile $outdir/landick_sig70_1398_SE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1400_Sig70.sort.bam --outfile $outdir/landick_sig70_1400_SE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1402_Sig70.sort.bam --outfile $outdir/landick_sig70_1402_SE_SCC.txt --scc.range $shift --ncluster 12

# # landick rif sig70 - PE ChIP-seq
# dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_PET/rif_treatment

# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1396_Sig70.sort.bam --outfile $outdir/landick_sig70_1396_PE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1398_Sig70.sort.bam --outfile $outdir/landick_sig70_1398_PE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1400_Sig70.sort.bam --outfile $outdir/landick_sig70_1400_PE_SCC.txt --scc.range $shift --ncluster 12
# rscripts/scripts/spp_QC_run.R --chipfile $dr/edsn1402_Sig70.sort.bam --outfile $outdir/landick_sig70_1402_PE_SCC.txt --scc.range $shift --ncluster 12





# # # landick aero sig70 - SE ChIP-seq
# # dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET/aerobic_vs_anaerobic_bed
# # isPET=FALSE

# # rscripts/scripts/SCC_curve_and_QC_ChIPseq_SE.R $dr/run80.sigma70+O2_I_PA.s_6_sequence.eland.Umatch.txt.bed $outdir/landick_sig70+O2_SE_SCC.txt $outdir2/landick_sig70+O2_SE_QC.txt $sizefile $isPET $maxShift

# # rscripts/scripts/SCC_curve_and_QC_ChIPseq_SE.R $dr/run80.sigma70-O2_IP_A.s_6_sequence.txt.eland.Umatch.txt.bed $outdir/landick_sig70-O2_SE_SCC.txt $outdir2/landick_sig70-O2_SE_QC.txt $sizefile $isPET $maxShift


# # # landick aero sig70 - PE ChIP-seq
# # dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_PET/aerobic_vs_anaerobic_bed
# # isPET=TRUE

# # rscripts/scripts/spp_QC_run.R $dr/run208.sigma70-minuso2a_TGACCA_L003.paired_Umatch_eland.tab.bed $dr/run208.sigma70-minuso2a_TGACCA_L003.paired_Umatch_eland.tab_coord.txt $outdir/landick_sig70-O2_PE_SCC.txt $outdir2/landick_sig70-O2_PE_QC.txt $sizefile $isPET $maxShift

# # rscripts/scripts/C_ChIPseq_PE.R $dr/run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab.bed $dr/run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab_coord.txt $outdir/landick_sig70+O2_PE_SCC.txt $outdir2/landick_sig70+O2_PE_QC.txt $sizefile $isPET $maxShift


