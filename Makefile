
# Loads the bam files into RData file format
load:
	R CMD BATCH inst/scripts/Script_loadData.R &

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

clean_out:
	rm -fr *Rout

# documents some functions
doc:
	R -e 'roxygen2::roxygenize()'

# conditional density
cond_density:
	R CMD BATCH inst/scripts/Script_conditionalDensity.R 

# pbc bottleneck coeff
pbc:
	R CMD BATCH inst/scripts/Script_PCR_bottleneck_coeff.R

# knit the vignettes
vignettes/%.md:vignettes/%.Rmd
	cd vignettes;R -e 'library(knitr);knit("$(<F)")';cd ..
