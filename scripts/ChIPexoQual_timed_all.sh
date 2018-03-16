#!/bin/sh

# This script


# Base

cores=4
outdir=./data/timed_runs/cores"$cores"

basedir=/p/keles/ChIPexo/volume4

CTCFdir=$basedir/pugh_data/mace
ERdir=$basedir/carroll_data/human/bamfiles
GRdir=$basedir/meijsing_data
TBPdir1=$basedir/venters_data/sortbam
TBPdir2=$basedir/zeitlinger_data/bam/sortbam


Rcall=./rscripts/scripts/ChIPexoQual_timed_run.R

# CTCF
$Rcall --bamfile $CTCFdir/CTCF_replicate1.sorted.bam --timefile $outdir/CTCF_rep1_time.tsv --mc.cores $cores
$Rcall --bamfile $CTCFdir/CTCF_replicate1.sorted.bam --timefile $outdir/CTCF_rep2_time.tsv --mc.cores $cores

# ER
$Rcall --bamfile $ERdir/ERR336950.sort.bam --timefile $outdir/ER_rep2_time.tsv --mc.cores $cores

# GR
$Rcall --bamfile $GRdir/IMR90_GR_chip-exo.sort.bam --timefile $outdir/GR_IMR90_time.tsv --mc.cores $cores

# TBPexo
$Rcall --bamfile $TBPdir1/TBP_K562_Rep3.sort.bam --timefile $outdir/TBPexo_K562_time.tsv --mc.cores $cores

# TBPnexus
$Rcall --bamfile $TBPdir2/ChIPnexus_K562_TBP_rep2.sort.bam --timefile $outdir/TBPnexus_K562_time.tsv --mc.cores $cores
