#!/bin/sh

outdir=data/SCC_spp_curves
outdir2=data/QC_spp
maxShift=300
cores=12

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336942.sort.bam --sccfile $outdir/carroll_mouse_Rep1_SCC.txt --outfile $outdir2/carroll_mouse_Rep1_SPP.txt --ncluster 12 
rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336956.sort.bam --sccfile $outdir/carroll_mouse_Rep2_SCC.txt --outfile $outdir2/carroll_mouse_Rep2_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336935.sort.bam --sccfile $outdir/carroll_mouse_Rep3_SCC.txt --outfile $outdir2/carroll_mouse_Rep3_SPP.txt --ncluster 12

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336933.sort.bam --sccfile $outdir/carroll_human_Rep1_SCC.txt --outfile $outdir2/carroll_human_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336950.sort.bam --sccfile $outdir/carroll_human_Rep2_SCC.txt --outfile $outdir2/carroll_human_Rep2_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ERR336938.sort.bam --sccfile $outdir/carroll_human_Rep3_SCC.txt --outfile $outdir2/carroll_human_Rep3_SPP.txt --ncluster 12

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data

rscripts/scripts/SPP_PBC.R --chipfile $dr/CTCF.bam --sccfile $outdir/pugh_CTCF_SCC.txt --outfile $outdir2/pugh_CTCF_SPP.txt --ncluster 12

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam --sccfile $outdir/chipnexus_embryo_dorsal_Rep1_SCC.txt --outfile $outdir2/chipnexus_embryo_dorsal_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam --sccfile $outdir/chipnexus_embryo_dorsal_Rep2_SCC.txt --outfile $outdir2/chipnexus_embryo_dorsal_Rep2_SPP.txt--ncluster 12 

rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_embryo_twist_rep1.sort.bam --sccfile $outdir/chipnexus_embryo_twist_Rep1_SCC.txt --outfile $outdir2/chipnexus_embryo_twist_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_embryo_twist_rep2.sort.bam --sccfile $outdir/chipnexus_embryo_twist_Rep2_SCC.txt --outfile $outdir2/chipnexus_embryo_twist_Rep2_SPP.txt --ncluster 12


rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_K562_TBP_rep1.sort.bam --sccfile $outdir/chipnexus_K562_TBP_Rep1_SCC.txt --outfile $outdir2/chipnexus_K562_TBP_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_K562_TBP_rep2.sort.bam --sccfile $outdir/chipnexus_K562_TBP_Rep2_SCC.txt --outfile $outdir2/chipnexus_K562_TBP_Rep2_SPP.txt --ncluster 12


rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_S2_Max_rep1.sort.bam --sccfile $outdir/chipnexus_S2_Max_Rep1_SCC.txt --outfile $outdir2/chipnexus_S2_Max_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_S2_Max_rep2.sort.bam --sccfile $outdir/chipnexus_S2_Max_Rep2_SCC.txt --outfile $outdir2/chipnexus_S2_Max_Rep2_SPP.txt --ncluster 12


rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_S2_MyC_rep1.sort.bam --sccfile $outdir/chipnexus_S2_MyC_Rep1_SCC.txt --outfile $outdir2/chipnexus_S2_MyC_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/ChIPnexus_S2_MyC_rep2.sort.bam --sccfile $outdir/chipnexus_S2_MyC_Rep2_SCC.txt --outfile $outdir2/chipnexus_S2_MyC_Rep2_SPP.txt --ncluster 12


# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/SPP_PBC.R --chipfile $dr/TBP_K562_Rep1.sort.bam --sccfile $outdir/venters_TBP_K562_Rep1_SCC.txt --outfile $outdir2/venters_TBP_K562_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/TBP_K562_Rep2.sort.bam --sccfile $outdir/venters_TBP_K562_Rep2_SCC.txt --outfile $outdir2/venters_TBP_K562_Rep2_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/TBP_K562_Rep3.sort.bam --sccfile $outdir/venters_TBP_K562_Rep3_SCC.txt --outfile $outdir2/venters_TBP_K562_Rep3_SPP.txt --ncluster 12

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

rscripts/scripts/SPP_PBC.R --chipfile $dr/IMR90_GR_chip-exo.sort.bam --sccfile $outdir/meijsing_GR_IMR90_Rep1_SCC.txt --outfile $outdir2/meijsing_GR_IMR90_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/K562_GR_chip-exo.sort.bam --sccfile $outdir/meijsing_GR_K562_Rep1_SCC.txt --outfile $outdir2/meijsing_GR_K562_Rep1_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/U2OS_GR_chip-exo.sort.bam --sccfile $outdir/meijsing_GR_U2OS_Rep1_SCC.txt --outfile $outdir2/meijsing_GR_U2OS_Rep1_SPP.txt --ncluster 12


# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic

rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn931_Sig70.sort.bam --sccfile $outdir/landick_sig70_931_SCC.txt --outfile $outdir2/landick_sig70_931_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn933_Sig70.sort.bam --sccfile $outdir/landick_sig70_933_SCC.txt --outfile $outdir2/landick_sig70_933_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn935_Sig70.sort.bam --sccfile $outdir/landick_sig70_935_SCC.txt --outfile $outdir2/landick_sig70_935_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn937_Sig70.sort.bam --sccfile $outdir/landick_sig70_937_SCC.txt --outfile $outdir2/landick_sig70_937_SPP.txt --ncluster 12

# landick aerobic sig70 rif
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn1311_Sig70.sort.bam --sccfile $outdir/landick_sig70_1311_SCC.txt --outfile $outdir2/landick_sig70_1311_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn1314_Sig70.sort.bam --sccfile $outdir/landick_sig70_1314_SCC.txt --outfile $outdir2/landick_sig70_1314_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn1317_Sig70.sort.bam --sccfile $outdir/landick_sig70_1317_SCC.txt --outfile $outdir2/landick_sig70_1317_SPP.txt --ncluster 12
rscripts/scripts/SPP_PBC.R --chipfile $dr/edsn1320_Sig70.sort.bam --sccfile $outdir/landick_sig70_1320_SCC.txt --outfile $outdir2/landick_sig70_1320_SPP.txt --ncluster 12
