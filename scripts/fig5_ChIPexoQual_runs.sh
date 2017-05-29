#!/bin/sh

pipeline=./rscripts/figures/fig5/fig5_ChIPexoQual_pipeline.R
outdr=./data/figures/fig5
cores=22
nregions=1000
ntimes=1000
basedr=/p/keles/ChIPexo/volume4

#########################################################

# FoxA1 mouse - carroll
indr=$basedr/carroll_data/mouse

$pipeline --readsfile $indr/ERR336935.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_FoxA1_mouseliver_Rep3_stats.tsv --scoresfile $outdr/carroll_FoxA1_mouseliver_Rep3_scores.tsv --cores $cores
$pipeline --readsfile $indr/ERR336942.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_FoxA1_mouseliver_Rep1_stats.tsv --scoresfile $outdr/carroll_FoxA1_mouseliver_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ERR336956.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_FoxA1_mouseliver_Rep2_stats.tsv --scoresfile $outdr/carroll_FoxA1_mouseliver_Rep2_scores.tsv --cores $cores 


#########################################################

# ER human - carroll

indr=$basedr/carroll_data/human/bamfiles

$pipeline --readsfile $indr/ERR336933.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_ER_MCF7_Rep1_stats.tsv --scoresfile $outdr/carroll_ER_MCF7_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ERR336950.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_ER_MCF7_Rep2_stats.tsv --scoresfile $outdr/carroll_ER_MCF7_Rep2_scores.tsv --cores $cores
$pipeline --readsfile $indr/ERR336938.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/carroll_ER_MCF7_Rep3_stats.tsv --scoresfile $outdr/carroll_ER_MCF7_Rep3_scores.tsv --cores $cores

#########################################################


# pugh original ctcf
indr=$basedr/pugh_data/mace

$pipeline --readsfile $indr/CTCF_replicate1.sorted.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_CTCF_HeLa_Rep1_stats.tsv --scoresfile $outdr/pugh_CTCF_HeLa_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/CTCF_replicate2.sorted.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_CTCF_HeLa_Rep2_stats.tsv --scoresfile $outdr/pugh_CTCF_HeLa_Rep2_scores.tsv --cores $cores
$pipeline --readsfile $indr/CTCF_replicate3.sorted.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_CTCF_HeLa_Rep3_stats.tsv --scoresfile $outdr/pugh_CTCF_HeLa_Rep3_scores.tsv --cores $cores

#########################################################

# chip nexus zeitlinger
indr=$basedr/zeitlinger_data/bam/sortbam

$pipeline --readsfile $indr/ChIPnexus_embryo_dorsal_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_dorsal_embryo_Rep1_stats.tsv --scoresfile $outdr/ChIPnexus_dorsal_embryo_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ChIPnexus_embryo_dorsal_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_dorsal_embryo_Rep2_stats.tsv --scoresfile $outdr/ChIPnexus_dorsal_embryo_Rep2_scores.tsv --cores $cores

$pipeline --readsfile $indr/ChIPnexus_embryo_twist_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_twist_embryo_Rep1_stats.tsv --scoresfile $outdr/ChIPnexus_twist_embryo_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ChIPnexus_embryo_twist_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_twist_embryo_Rep2_stats.tsv --scoresfile $outdr/ChIPnexus_twist_embryo_Rep2_scores.tsv --cores $cores

$pipeline --readsfile $indr/ChIPnexus_K562_TBP_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_TBP_K562_Rep1_stats.tsv --scoresfile $outdr/ChIPnexus_TBP_K562_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ChIPnexus_K562_TBP_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_TBP_K562_Rep2_stats.tsv --scoresfile $outdr/ChIPnexus_TBP_K562_Rep2_scores.tsv --cores $cores

$pipeline --readsfile $indr/ChIPnexus_S2_Max_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_Max_S2_Rep1_stats.tsv --scoresfile $outdr/ChIPnexus_Max_S2_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ChIPnexus_S2_Max_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_Max_S2_Rep2_stats.tsv --scoresfile $outdr/ChIPnexus_Max_S2_Rep2_scores.tsv --cores $cores

$pipeline --readsfile $indr/ChIPnexus_S2_MyC_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_MyC_S2_Rep1_stats.tsv --scoresfile $outdr/ChIPnexus_MyC_S2_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/ChIPnexus_S2_MyC_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/ChIPnexus_MyC_S2_Rep2_stats.tsv --scoresfile $outdr/ChIPnexus_MyC_S2_Rep2_scores.tsv --cores $cores

#########################################################

# venters TBP
indr=$basedr/venters_data/sortbam

$pipeline --readsfile $indr/TBP_K562_Rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/venters_TBP_K562_Rep1_stats.tsv --scoresfile $outdr/venters_TBP_K562_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/TBP_K562_Rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/venters_TBP_K562_Rep2_stats.tsv --scoresfile $outdr/venters_TBP_K562_Rep2_scores.tsv --cores $cores
$pipeline --readsfile $indr/TBP_K562_Rep3.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/venters_TBP_K562_Rep3_stats.tsv --scoresfile $outdr/venters_TBP_K562_Rep3_scores.tsv --cores $cores

#########################################################

# meijsing GR
indr=$basedr/meijsing_data/

$pipeline --readsfile $indr/IMR90_GR_chip-exo.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/meijsing_GR_IMR90_Rep1_stats.tsv --scoresfile $outdr/meijsing_GR_IMR90_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/K562_GR_chip-exo.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/meijsing_GR_K562_Rep1_stats.tsv --scoresfile $outdr/meijsing_GR_K562_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/U2OS_GR_chip-exo.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/meijsing_GR_U2OS_Rep1_stats.tsv --scoresfile $outdr/meijsing_GR_U2OS_Rep1_scores.tsv --cores $cores

#########################################################

# exo histone

indr=$basedr/exo_histone_data

$pipeline --readsfile $indr/BAM/histone1_rep1.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_histone1_yeast_Rep1_stats.tsv --scoresfile $outdr/pugh_histone1_yeast_Rep1_scores.tsv --cores $cores
$pipeline --readsfile $indr/BAM/histone1_rep2.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_histone1_yeast_Rep2_stats.tsv --scoresfile $outdr/pugh_histone1_yeast_Rep2_scores.tsv --cores $cores
$pipeline --readsfile $indr/BAM/histone1_rep3.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_histone1_yeast_Rep3_stats.tsv --scoresfile $outdr/pugh_histone1_yeast_Rep3_scores.tsv --cores $cores
$pipeline --readsfile $indr/BAM/histone1_rep4.sort.bam --nregions $nregions --ntimes $ntimes --statsfile $outdr/pugh_histone1_yeast_Rep4_stats.tsv --scoresfile $outdr/pugh_histone1_yeast_Rep4_scores.tsv --cores $cores
