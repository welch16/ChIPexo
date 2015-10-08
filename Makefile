
# Loads the bam files into RData file format
load:
	R CMD BATCH inst/rscripts/Script_loadData.R &

# removes not necessary stuff
clean:
	rm -fr *~
	rm -fr */*~
	rm -f .RDataTmp
	rm -f .Rhistory
	rm -f */.Rhistory
	rm -f inst/rscripts/*~
	rm -f inst/rscripts/Ren/*~
	rm -f inst/rscripts/*Rout
	rm -f .RData
	rm -f .*~

clean_out:
	rm -fr *Rout

# documents some functions
doc:
	R -e 'roxygen2::roxygenize()'

# conditional density
cond_density:
	R CMD BATCH inst/rscripts/Script_conditionalDensity.R 

# densities
densities:
	R CMD BATCH inst/rscripts/Script_density.R


# pbc bottleneck coeff
pbc:
	R CMD BATCH inst/rscripts/Script_PCR_bottleneck_coeff.R

# depth scripth and logcount boxplots
depth:
	R CMD BATCH inst/rscripts/Script_depth.R 

# cross-correlation 
cross:
	R CMD BATCH inst/rscripts/Script_crossCorr_exploration.R

# peak-plots
peaks_plots:
	R CMD BATCH inst/rscripts/Script_cluster_structure.R 

peaks_slice:
	R CMD BATCH inst/rscripts/Script_cluster_structure_slice.R 

peaks_slice_Ren:
	R CMD BATCH inst/rscripts/Ren/Script_cluster_analysis.R

# diff-plots
diff:
	R CMD BATCH inst/rscripts/Script_detailled_summary.R


# positions
positions:
	R CMD BATCH inst/rscripts/Script_positions.R

# knit the vignettes
vignettes/%.md:vignettes/%.Rmd
	cd vignettes;R -e 'library(knitr);knit("$(<F)")';cd ..

vignettest/%.pdf:vignettes/%.tex
	cd vignettes;pdflatex $(<F);bibtex $(<F);pdflatex $(<F);pdflatex $(<F);cd ..

vignettes/%.tex:vignettes/%.Rnw
	cd vignettes;R CMD Sweave --engine=knitr::knitr --pdf $(<F);cd ..

sample_saturation:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/sample_chip.R

bins:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_exo_construct_bins.R

bins_pet:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_seq_pet_construct_bins.R

bins_set:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_seq_set_construct_bins.R

peaks:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_exo_call_peaks.R

peaks_pet:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_seq_pet_call_peaks.R

peaks_set:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_seq_set_call_peaks.R

binding_sites:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/chip_exo_call_bs.R

sample_SET:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save rscripts/create_ChIPseq_SET.R
