
# Loads the bam files into RData file format
load:
	R CMD BATCH R/Script_loadData.R &

clean:
	rm -fr *~
	rm -fr */*~
	rm -f .RDataTmp
	rm -f .Rhistory
	rm -f */.Rhistory
