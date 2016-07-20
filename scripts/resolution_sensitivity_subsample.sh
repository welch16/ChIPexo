#!/bin/sh

indir=/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity

nreads=7400000
npairs=`expr $nreads / 2`
seed=12345

# generate indexes

here=$PWD

cd $indir

## ChIP-exo rep 1
cd ChIPexo

perl ../index_reads_SET.pl ChIP-exo_sigma70_exp_phase_R1.sam index_1_$seed.txt $nreads 2 $seed
wait
perl ../subsample_reads.pl ChIP-exo_sigma70_exp_phase_R1.sam index_1_$seed.txt & 

## ChIP-exo rep 2

perl ../index_reads_SET.pl ChIP-exo_sigma70_exp_phase_R2.sam index_2_$seed.txt $nreads 2 $seed
wait
perl ../subsample_reads.pl ChIP-exo_sigma70_exp_phase_R2.sam index_2_$seed.txt &

cd ..

## ChIP-seq PET

cd ChIPseq_PET

perl ../index_reads_PET.pl run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab index_pet_$seed.txt $npairs 2 $seed
wait
perl ../subsample_reads.pl run208.sigma70-pluso2a_TTAGGC_L003.paired_Umatch_eland.tab index_pet_$seed.txt &

cd ..

## ChIP-seq SET

cd ChIPseq_SET

perl ../index_reads_SET.pl run80.sigma70+O2_I_PA.s_6_sequence.eland.Umatch.txt index_set_$seed.txt $nreads 2 $seed 
wait
perl ../subsample_reads.pl run80.sigma70+O2_I_PA.s_6_sequence.eland.Umatch.txt index_set_$seed.txt & 

cd $here
