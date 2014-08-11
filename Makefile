
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

# documents some functions
doc:
	R -e 'roxygen2::roxygenize()'
