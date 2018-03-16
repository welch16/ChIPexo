
## ChIP-exo , PE ChIP-seq and SE ChIP-seq comparison

Wanna compare in resolution and sensitivity by sampling 6M reads from
ChIP-exo and SE ChIP-seq. For PE ChIP-Seq want to sample 3M pairs
(i.e. 6M reads).

### Resolution and sensitivity

[Code](/p/keles/ChIPexo/volume3/ChIPexo/paper_analysis/resolution_sensitivity_R1_old_data.R)

Requires:

		* Peaks
		* Binding sites estimation
		* Promoters table

### Promoters table

[Promoters](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/database_regulonDB/PromoterSigma70Set.txt)

### Binding sites estimation

[ChIP-exo Rep1](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R1/dpeak_ChIP-exo_R1.R)
[ChIP-exo Rep2](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R1/dpeak_ChIP-exo_R1.R)

[ChIP-seq PE](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_PET/sigma70+O2/mosaics_PET/sigma70+O2_decon_PET_Kiley.R)

[ChIP-seq SE](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_SET/sigma70+O2/mosaics_SET/sigma70+O2_analysis.R)

Requires:

		* Peaks 
		* Reads

### Peaks

[ChIP-exo](/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_mosaics/mosaics_ChIP-exo_1S_new_R1.R)

[ChIP-seq PE]()
[ChIP-seq SE]()
	
Requires:

		* BinData -> Reads
		* GC, M, N
		* Input samples



