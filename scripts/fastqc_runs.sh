#!/bin/sh

outdr=data/fastqc_runs

# carroll mouse
dr=/p/keles/ChIPexo/volume4/carroll_data/mouse

/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336942.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336956.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336935.sort.bam &

# carroll human
dr=/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles

/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336933.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336950.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ERR336938.sort.bam & 

# pugh original ctcf
dr=/p/keles/ChIPexo/volume4/pugh_data/raw_data

/p/stat/genomics/bin/fastqc -o $outdr $dr/CTCF.bam 

# chip nexus zeitlinger
dr=/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam 

/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_embryo_dorsal_rep1.sort.bam & 
/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_embryo_dorsal_rep2.sort.bam &

/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_embryo_twist_rep1.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_embryo_twist_rep2.sort.bam &


/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_K562_TBP_rep1.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_K562_TBP_rep2.sort.bam &


/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_S2_Max_rep1.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_S2_Max_rep2.sort.bam &


/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_S2_MyC_rep1.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/ChIPnexus_S2_MyC_rep2.sort.bam 


# venters TBP
dr=/p/keles/ChIPexo/volume4/venters_data/sortbam

/p/stat/genomics/bin/fastqc -o $outdr $dr/TBP_K562_Rep1.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/TBP_K562_Rep2.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/TBP_K562_Rep3.sort.bam &

# meijsing GR
dr=/p/keles/ChIPexo/volume4/meijsing_data/

/p/stat/genomics/bin/fastqc -o $outdr $dr/IMR90_GR_chip-exo.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/K562_GR_chip-exo.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/U2OS_GR_chip-exo.sort.bam


# landick aerobic sig70
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic

/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn931_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn933_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn935_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn937_Sig70.sort.bam &

# landick aerobic sig70 rif
dr=/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment

/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn1311_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn1314_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn1317_Sig70.sort.bam &
/p/stat/genomics/bin/fastqc -o $outdr $dr/edsn1320_Sig70.sort.bam &

