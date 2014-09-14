
# Loads the bam files into RData file format
load:
	R CMD BATCH inst/rscripts/Script_loadData.R &
	R CMD BATCH inst/rscripts/Script_loadRenData.R &

# removes not necessary stuff
clean:
	rm -fr *~
	rm -fr */*~
	rm -f .RDataTmp
	rm -f .Rhistory
	rm -f */.Rhistory
	rm -f inst/scripts/*~
	rm -f inst/scripts/*Rout
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

# knit the vignettes
vignettes/%.md:vignettes/%.Rmd
	cd vignettes;R -e 'library(knitr);knit("$(<F)")';cd ..

