#!/bin/sh

outdir2=data/QC_htseqtools
cores=16

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336942.sort.bam --outfile $outdir2/carroll_mouse_Rep1_htseqtools.txt --ncluster $(cores) 
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336956.sort.bam --outfile $outdir2/carroll_mouse_Rep2_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336935.sort.bam --outfile $outdir2/carroll_mouse_Rep3_htseqtools.txt --ncluster $(cores)

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336933.sort.bam --outfile $outdir2/carroll_human_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336950.sort.bam --outfile $outdir2/carroll_human_Rep2_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ERR336938.sort.bam --outfile $outdir2/carroll_human_Rep3_htseqtools.txt --ncluster $(cores)

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data/raw_data

rscripts/scripts/htseqtools_qc.R --chipfile $dr/CTCF.bam --outfile $outdir2/pugh_CTCF_htseqtools.txt --ncluster $(cores)

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam

rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam --outfile $outdir2/chipnexus_embryo_dorsal_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam --outfile $outdir2/chipnexus_embryo_dorsal_Rep2_htseqtools.txt--ncluster $(cores)

rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_embryo_twist_rep1.sort.bam --outfile $outdir2/chipnexus_embryo_twist_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_embryo_twist_rep2.sort.bam --outfile $outdir2/chipnexus_embryo_twist_Rep2_htseqtools.txt --ncluster $(cores)


rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_K562_TBP_rep1.sort.bam --outfile $outdir2/chipnexus_K562_TBP_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_K562_TBP_rep2.sort.bam --outfile $outdir2/chipnexus_K562_TBP_Rep2_htseqtools.txt --ncluster $(cores)


rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_S2_Max_rep1.sort.bam --outfile $outdir2/chipnexus_S2_Max_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_S2_Max_rep2.sort.bam --outfile $outdir2/chipnexus_S2_Max_Rep2_htseqtools.txt --ncluster $(cores)


rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_S2_MyC_rep1.sort.bam --outfile $outdir2/chipnexus_S2_MyC_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/ChIPnexus_S2_MyC_rep2.sort.bam --outfile $outdir2/chipnexus_S2_MyC_Rep2_htseqtools.txt --ncluster $(cores)


# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

rscripts/scripts/htseqtools_qc.R --chipfile $dr/TBP_K562_Rep1.sort.bam --outfile $outdir2/venters_TBP_K562_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/TBP_K562_Rep2.sort.bam --outfile $outdir2/venters_TBP_K562_Rep2_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/TBP_K562_Rep3.sort.bam --outfile $outdir2/venters_TBP_K562_Rep3_htseqtools.txt --ncluster $(cores)

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

rscripts/scripts/htseqtools_qc.R --chipfile $dr/IMR90_GR_chip-exo.sort.bam --outfile $outdir2/meijsing_GR_IMR90_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/K562_GR_chip-exo.sort.bam --outfile $outdir2/meijsing_GR_K562_Rep1_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/U2OS_GR_chip-exo.sort.bam --outfile $outdir2/meijsing_GR_U2OS_Rep1_htseqtools.txt --ncluster $(cores)


# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic

rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn931_Sig70.sort.bam --outfile $outdir2/landick_sig70_931_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn933_Sig70.sort.bam --outfile $outdir2/landick_sig70_933_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn935_Sig70.sort.bam --outfile $outdir2/landick_sig70_935_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn937_Sig70.sort.bam --outfile $outdir2/landick_sig70_937_htseqtools.txt --ncluster $(cores)

# landick aerobic sig70 rif
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn1311_Sig70.sort.bam --outfile $outdir2/landick_sig70_1311_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn1314_Sig70.sort.bam --outfile $outdir2/landick_sig70_1314_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn1317_Sig70.sort.bam --outfile $outdir2/landick_sig70_1317_htseqtools.txt --ncluster $(cores)
rscripts/scripts/htseqtools_qc.R --chipfile $dr/edsn1320_Sig70.sort.bam --outfile $outdir2/landick_sig70_1320_htseqtools.txt --ncluster $(cores)
