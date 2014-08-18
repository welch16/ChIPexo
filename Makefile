
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
	rm -f *.Rout

# documents some functions
doc:
	R -e 'roxygen2::roxygenize()'

# conditional density
cond_density:
	R CMD BATCH inst/scripts/Script_conditionalDensity.R 
