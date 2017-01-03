#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

## defines option parse and help caption for script
optList = list(
  make_option("--chipfile" , action = "store_true",type = "character",
              help = "File corresponding to the ChIP-seq/exo/nexus file used"),
  make_option("--outfile", action = "store_true",type = "character",default = tempfile(),
              help = "File where the SCC is going to be stored"),
  make_option("--scc.range",action = "store_true",type = "character",default = "1,300",
              help = "Range used to calculate the SCC separated by a comma ','."),
  make_option("--ncluster",action = "store_true",type = "numeric",default = 12,
              help = "Number of clusters used for parallel processing.")
 )


opt = parse_args(OptionParser(option_list = optList))

## opt$chipfile = "/p/keles/ChIPexo/volume4/carroll_data/mouse/ERR336956.sort.bam"

stopifnot(file.exists(opt$chipfile))

library(spp)
library(snow)
library(parallel)
library(readr)

sccrange = as.numeric(strsplit(opt$scc.range,"\\,")[[1]])

tags = read.bam.tags(opt$chipfile)

cl = makeCluster(opt$ncluster, type = "SOCK")

scc = get.binding.characteristics(tags,srange = sccrange,bin = 1,cluster = cl)

cross = scc$cross.correlation
names(cross) = c("shift","cross.corr")
cross$shift = as.integer(cross$shift)

write_delim(cross,path = opt$outfile,delim = "\t")


