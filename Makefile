
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
peaks:
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

vignettes/%.pdf:vignettes/%.Rnw
	cd vignettes;R CMD Sweave --engine=knitr::knitr --pdf $(<F);cd ..
