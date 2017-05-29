#!/bin/sh

# This is a quick run of the FASTQC pipeline
# on the FoxA1 mouse liver samples

in_dir=/p/keles/ChIPexo/volume4/carroll_data/mouse
out_dir=otherQC_methods/FASTQC

rep3=ERR336935.sort.bam
rep1=ERR336942.sort.bam
rep2=ERR336956.sort.bam

# mkdir $out_dir/FoxA1_Rep1
# mkdir $out_dir/FoxA1_Rep2
# mkdir $out_dir/FoxA1_Rep3

# fastqc -o $out_dir/FoxA1_Rep1 -f bam -t 8 $in_dir/$rep1 &
# fastqc -o $out_dir/FoxA1_Rep2 -f bam -t 8 $in_dir/$rep2 & 
# fastqc -o $out_dir/FoxA1_Rep3 -f bam -t 8 $in_dir/$rep3 &


in_dir=/p/keles/ChIPexo/volume4/venters_data/sortbam

rep1=TBP_K562_Rep1.sort.bam
rep2=TBP_K562_Rep2.sort.bam
rep3=TBP_K562_Rep3.sort.bam

# mkdir $out_dir/TBP_K562_exo_rep1
# mkdir $out_dir/TBP_K562_exo_rep2
# mkdir $out_dir/TBP_K562_exo_rep3

# fastqc -o $out_dir/TBP_K562_exo_rep1 -f bam -t 8 $in_dir/$rep1 &
# fastqc -o $out_dir/TBP_K562_exo_rep2 -f bam -t 8 $in_dir/$rep2 &
# fastqc -o $out_dir/TBP_K562_exo_rep3 -f bam -t 8 $in_dir/$rep3 &

in_dir=/p/keles/ChIPexo/volume4/zeitlinger_data/bam

rep1=ChIPnexus_K562_TBP_rep1.bam
rep2=ChIPnexus_K562_TBP_rep2.bam

mkdir $out_dir/TBP_K562_nexus_rep1
mkdir $out_dir/TBP_K562_nexus_rep2

fastqc -o $out_dir/TBP_K562_nexus_rep1 -f bam -t 8 $in_dir/$rep1 & 
fastqc -o $out_dir/TBP_K562_nexus_rep2 -f bam -t 8 $in_dir/$rep2 & 
